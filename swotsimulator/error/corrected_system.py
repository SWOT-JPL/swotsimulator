# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Estimated System errors
---------------------
"""
from typing import Dict, Tuple
import numpy as np
import netCDF4

from .. import settings

class CorrectedSystem:
    """
    Corrected system errors

    Args:
        parameters (settings.Parameters): Simulation settings
        first_date (numpy.datetime64): Date of the first simulated
            measurement.
    """
    def __init__(self, parameters: settings.Parameters,
                 first_date: np.datetime64,
                 list_error: Optional[list] = None) -> None:
        self.pattern = parameters.corrected_system_dataset
        if list_error is None:
            self.listerr = parameters.list_error
        else:
            self.listerr = list_error
        self.proll_err = ds['proll_err'][:]
        self.p1phase_err = ds['p1phase_err'][:]
        self.p2phase_err = ds['p2phase_err'][:]
        self.slope1_est = ds['slope1_est'][:]
        self.slope2_est = ds['slope2_est'][:]
        self.slope1_err = ds['slope1_err'][:]
        self.slope2_err = ds['slope2_err'][:]
        ds.close()

    def _get_file_name(self, ipas: int, cycle: int) -> str:
        ifile = f'{self.pattern}_c{cycle:03d}_p{ipass:04d}.nc'
        return ifile

    def _read_data(self, ifile: str, listerr: list
                  ) -> Tuple[np.ndarray, dict]:
        ds = netCDF4.Dataset(ifile, 'r')
        _time_date = first_date + (ds['time'][:].astype(np.float32)
                     * 1000000).astype("timedelta64[us]")
        time_date = (time_date.astype('datetime64[us]').astype('float64')
                          * 0.001)
        dic_err = {}
        for ivar in self.listerr:
            dic_err[ivar] = ds[ivar][:]
        ds.close()
        return time_date, dic_err

    def _generate_1d(self, time: np.ndarray
                     ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        ifile = self._get_file_name()
        _time_date, _system_err = self._read_data(ifile)
        time = time.astype('datetime64[us]').astype('float64')
        phase_l = np.interp(time, self.time_date, self.p1phase_err)
        phase_r = np.interp(time, self.time_date, self.p2phase_err)
        phase1d = np.concatenate(([phase_l.T], [phase_r.T]), axis=0)
        roll = np.interp(time, self.time_date, self.proll_err)

        est_l = np.interp(time, self.time_date, self.slope1_est)
        err_l = np.interp(time, self.time_date, self.slope1_err)
        rem_l = est_l - err_l
        est_r = np.interp(time, self.time_date, self.slope2_est)
        err_r = np.interp(time, self.time_date, self.slope2_err)
        rem_r = est_r - err_r
        theta2 = np.concatenate(([rem_l.T], [rem_r.T]), axis=0)

        return roll * 1e-3, phase1d.T * 1e-3, theta2.T * 1e-3

    def generate(
        self,
        time: np.ndarray,
        x_ac: np.ndarray,
    ) -> Dict[str, np.ndarray]:
        """Interpolate roll and phase and errors

        Args:
            time (numpy.ndarray): Date of measurements.
            x_ac (numpy.ndarray): Across track distance.

        Returns:
            dict: variable name and errors simulated.
        """
        roll_1d, phase_1d, rollphase_est_1d = self._generate_1d(time)
        num_pixels = x_ac.shape[0]
        swath_center = num_pixels // 2
        ac_l = x_ac[:swath_center]
        ac_r = x_ac[swath_center:]

        phase = np.full((phase_1d.shape[0], num_pixels), np.nan)
        phase[:, :swath_center] = ac_l * phase_1d[:, 0, np.newaxis]
        phase[:, swath_center:] = ac_r * phase_1d[:, 1, np.newaxis]

        rollphase_est = np.full((phase_1d.shape[0], num_pixels), np.nan)
        rollphase_est[:, :swath_center] = np.mat(rollphase_est_1d[:,
                                                                  0]).T * ac_l
        rollphase_est[:,
                      swath_center:] = np.mat(rollphase_est_1d[:, 1]).T * ac_r
        return {
            "simulated_error_roll": x_ac * roll_1d[:, np.newaxis],
            "simulated_error_phase": phase,
            "simulated_roll_phase_estimate": rollphase_est
        }
