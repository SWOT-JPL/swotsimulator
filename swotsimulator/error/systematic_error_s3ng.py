# Copyright (c) 2024 OceanDataLab
#
"""
Estimated Systematic errors for S3NG
------------------------------------
"""
from typing import Optional
import numpy as np
import netCDF4
from .. import settings
from scipy.interpolate import interp1d


class SystematicErrors3ng:
    """
    Systematic errors
    Args
    """

    def __init__(self, parameters: settings.Parameters,
                 first_date: np.datetime64,
                 list_error: Optional[list] = None) -> None:
        if list_error is None:
             self.listerror = ['all',]
        else:
             self.listerror = list_error
        first_time = (first_date
                 - np.datetime64( "2000-01-01 00:00:00")).astype("float")
        self.generate1d(parameters.file_systematic, first_time)

    def _read_data(self, filenc: str, first_time: float) -> None:

        with netCDF4.Dataset(filenc, 'r') as fcid:
            self.time_syst = first_time + np.array(fcid.variables['Time_second'][:])
            self.roll_gse = np.array(fcid.variables['baseline_gse_roll_microrad'][:])
            self.roll_ted = np.array(fcid.variables['baseline_ted_roll_microrad'][:])
            self.bd = np.array(fcid.variables['baseline_dilatation_micrometer'][:])
            self.phase_common = np.array(fcid.variables['common_diff_phase_degree'][:])
            self.phase_diff = np.array(fcid.variables['LR_diff_phase_degree'][:])

    def _convert(self) -> None:
        from .. import const_s3ng as const
        # Convert phase in microrad
        self.ephase = (const.C / (const.Fka * const.B)
                       * (1 + const.sat_elev / const.Rearth)
                       * (self.phase_common * np.pi / 180) * 1e6)  # in microrad
        self.ephased = (const.C / (const.Fka * const.B)
                        * (1 + const.sat_elev / const.Rearth)
                        * (self.phase_diff * np.pi / 180) * 1e6)  # in microrad

        # Convert baseline dilation in m-1
        self.hbd = ((1 + const.sat_elev / const.Rearth)
                    * 1 / (const.sat_elev * const.B)
                    * self.bd * 1e-6)  # in m / m**2

    def _interpolator(self, kind: Optional[str] = 'linear'):
        self.finterp_roll_gse = interp1d(self.time_syst, self.roll_gse,
                                         kind=kind, bounds_error=False)
        self.finterp_roll_ted = interp1d(self.time_syst, self.roll_ted,
                                         kind=kind, bounds_error=False)
        self.finterp_ephase = interp1d(self.time_syst, self.ephase,
                                       kind=kind, bounds_error=False)
        self.finterp_ephased = interp1d(self.time_syst, self.ephased,
                                        kind=kind, bounds_error=False)
        self.finterp_hbd = interp1d(self.time_syst, self.hbd,
                                    kind=kind, bounds_error=False)

    def generate1d(self, filenc: str, first_time: float) -> None:
        self._read_data(filenc, first_time)
        self._convert()
        self._interpolator()

    def generate(self, time: np.ndarray, x_ac: np.ndarray):
        num_pixels = x_ac.shape[0]
        # swath_center = int(num_pixels / 2)
        # ac_l = x_ac[:swath_center]
        # ac_r = x_ac[swath_center:]
        ntime = time.shape[0]
        x_acm = x_ac * 10**3
        tmp = self.finterp_ephase(time)
        ephase = np.full((ntime, num_pixels), np.nan)
        ephase[:, :] = 1e-6 * x_acm * tmp[:, np.newaxis]
        tmp = self.finterp_ephased(time)
        ephased = np.full((ntime, num_pixels), np.nan)
        ephased[:, :] = 1e-6 * abs(x_acm) * tmp[:, np.newaxis]
        tmp = self.finterp_roll_gse(time)
        roll_gse = np.full((ntime, num_pixels), np.nan)
        roll_gse = np.asmatrix(tmp).T * x_acm * 1e-6
        tmp = self.finterp_roll_ted(time)
        roll_ted = np.full((ntime, num_pixels), np.nan)
        roll_ted = np.asmatrix(tmp).T * x_acm * 1e-6
        tmp = self.finterp_hbd(time)
        hbd = np.full((ntime, num_pixels), np.nan)
        hbd = x_acm**2 * tmp[:, np.newaxis]
        return {"roll_gse": roll_gse,
                "roll_ted": roll_ted,
                "phase_relative": ephased,
                "phase_absolute": ephase,
                "baseline_dilation": hbd}
