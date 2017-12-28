'''Constants that are defined for SWOT simulator. \n
It contains Earth parameters as well as SWOT instrument
and satellite caracteristics. '''

# ################################
# # EARTH CONSTANTS             ##
# ################################
# - Earth radius (m)
Rearth = 6378. * 10**3
# - Convert degree to km
deg2km = 111.11
# - Seconds in a day
secinday = 86400.
# - Light speed (m/s)
C = 2.998*10**8


# ###################################
# # SWOT INSTRUMENT CARACTERISTICS ##
# ###################################
# - Satellite elevation (m)
sat_elev = 891*10**3
# - Baseline (m)
B = 10
# (in Hz)
Fka = 35.75 * 10**9
# Number of days of one SWOT cycle
tcycle = 20.86455

# ###################################
# # OTHER PARAMETERS               ##
# ###################################
# - Radius to interpolate locally model data on the swath (in km)
#  data are selected every xal_step points and on a radius of radius_interp
radius_interp = 100.
# - Sampling to interpolate locally model data on the swath (in km)
#  Data are selected every xal_step points and on a radius of radius_interp
xal_step = 20.

# ###################################
# # VARIABLE NAMES AND UNITS       ##
# ###################################
f8_nan = 1.36e9
i4_nan = 2147483647
u2_nan = 65536
u1_nan = 255
i2_nan = 32767
pd = {'varname': 'pd',
      'longname': 'Simulated path delay error due to wet tropo',
      'unit': 'm',
      'fill_value': i4_nan,
      'min_value': -1500000,
      'max_value': 1500000,
      'type': 'i4',
      'scale': 1e-4
      }
SSH_model = {'varname': 'ssh_model',
             'longname': 'SSH interpolated from model',
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }
SSH_obs = {'varname': 'ssh_obs',
           'longname': 'Observed SSH (ssh_model + errors)',
           'unit': 'm',
           'fill_value': i4_nan,
           'min_value': -1500000,
           'max_value': 1500000,
           'type': 'i4',
           'scale': 1e-4
           }
index = {'varname': 'model_index',
         'longname': 'Equivalent model output number in list of file',
         'unit': '',
         'min_value': 0,
         'max_value': 65535,
         'fill_value': u2_nan,
         'type': 'u2',
         'scale': 1
         }
longname = 'Residual path delay error after a 1-beam radiometer correction'
pd_err_1b = {'varname': 'pd_err_1b',
             'longname': longname,
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }

longname = 'Residual path delay error after a 2-beam radiometer correction'
pd_err_2b = {'varname': 'pd_err_2b',
             'longname': longname,
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }
roll_err = {'varname': 'roll_err',
            'longname': 'Residual roll error',
            'unit': 'm',
            'fill_value': i4_nan,
            'min_value': -1500000,
            'max_value': 1500000,
            'type': 'i4',
            'scale': 1e-4
            }
phase_err = {'varname': 'phase_err',
             'longname': 'Phase error',
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }
timing_err = {'varname': 'timing_err',
              'longname': 'Timing error',
              'unit': 'm',
              'fill_value': i4_nan,
              'min_value': -1500000,
              'max_value': 1500000,
              'type': 'i4',
              'scale': 1e-4
              }
bd_err = {'varname': 'bd_err',
          'longname': 'Baseline dilation error',
          'unit': 'm',
          'fill_value': i4_nan,
          'min_value': -1500000,
          'max_value': 1500000,
          'type': 'i4',
          'scale': 1e-4
          }
karin_err = {'varname': 'karin_err',
             'longname': 'Karin instrument random error',
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }
nadir_err = {'varname': 'nadir_err',
             'longname': 'Nadir altimeter error',
             'unit': 'm',
             'fill_value': i4_nan,
             'min_value': -1500000,
             'max_value': 1500000,
             'type': 'i4',
             'scale': 1e-4
             }

# - Variables for CNES mockup
list_var_mockup = ['ssha_karin_swath', 'ssha_uncert', 'ssh_quality_flag',
                   'karin_surf_type', 'dynamic_ice_flag', 'rad_surf_type',
                   'rain_flag', 'mss_reference', 'geoid',
                   'internal_tide_solution1', 'karin_karin_xover_height_corr',
                   'karin_na_height_corr', 'karin_cal_flag']

# Karin variables
longname = ('SSHA = SSH – MSS – ocean tide – load tide – solid tide – pole'
            'tide – inv_bar_corr – DAC_corr')
ssha_karin_swath = {'varname': 'SSHA_kaRIn_swath',
                    'longname': longname,
                    'unit': 'm',
                    'fill_value': i4_nan,
                    'min_value': -100000,
                    'max_value': 100000,
                    'type': 'i4',
                    'scale': 1e-4
                    }
ssha_uncert = {'varname': 'SSHA_uncert',
               'longname': 'Total SSHA uncertainty',
               'unit': 'm',
               'fill_value': u2_nan,
               'min_value': 0,
               'max_value': None,
               'type': 'u2',
               'scale': 1.e-4
               }
longname = 'Quality indicator on delivered SSH in the swath'
ssh_quality_flag = {'varname': 'SSH_quality_flag',
                    'longname': longname,
                    'unit': 'bits',
                    'fill_value': u1_nan,
                    'min_value': 0,
                    'max_value': 255,
                    'type': 'u1',
                    'scale': 1
                    }

# Geophysical flags
longname = ('KaRin surface type (0 : ocean, 1 : lakes/enclosed sea, 2 : ice, '
            '3 : land)')
karin_surf_type = {'varname': 'KaRIn_surf_type',
                   'longname': longname,
                   'unit': '',
                   'fill_value': u1_nan,
                   'min_value': 0,
                   'max_value': 3,
                   'type': 'u1',
                   'scale': 1,
                   }
dynamic_ice_flag = {'varname': 'dynamic_ice_flag',
                    'longname': 'Ice flag from EUMETSat and/or KaRIn data',
                    'unit': '',
                    'fill_value': u1_nan,
                    'min_value': 0,
                    'max_value': 1,
                    'type': 'u1',
                    'scale': 1,
                    }
longname = ('Radiometer surface type (0 : ocean, 1 : lakes/enclosed sea,'
            ' 2 : ice, 3 : land)')
rad_surf_type = {'varname': 'rad_surf_type',
                 'longname': longname,
                 'unit': '',
                 'fill_value': u1_nan,
                 'min_value': 0,
                 'max_value': 3,
                 'type': 'u1',
                 'scale': 1,
                 }
rain_flag = {'varname': 'rain_flag',
             'longname': 'Rain flag',
             'unit': 'bits',
             'fill_value': u1_nan,
             'min_value': 0,
             'max_value': 100,
             'type': 'u1',
             'scale': 1,
             }

# Geophysical references
longname = ('mean sea surface height above reference ellipsoid - Reference '
            'solution')
mss_reference = {'varname': 'mss_reference',
                 'longname': longname,
                 'unit': 'm',
                 'fill_value': i4_nan,
                 'min_value': -1500000,
                 'max_value': 1500000,
                 'type': 'i4',
                 'scale': 1e-4
                 }
geoid = {'varname': 'geoid',
	 'longname': 'geoid height above the ellipsoid',
	 'unit': 'm',
	 'fill_value': i4_nan,
	 'min_value': -1500000,
	 'max_value': 1500000,
	 'type': 'i4',
	 'scale': 1e-4
	 }
longname = 'Coherent internal tide model (Not applied)'
internal_tide_solution1 = {'varname': 'internal_tide_solution1',
                           'longname': longname,
			   'unit': 'm',
			   'fill_value': i2_nan,
			   'min_value': -2000,
			   'max_value': 2000,
			   'type': 'i2',
			   'scale': 1e-4
			   }
longname = 'Correction processed from KaRIn/KaRIn cross-overs'
karin_karin_xover_height_corr = {'varname': 'KaRIn_KaRIn_Xover_height_corr',
				 'longname': longname,
				 'unit': 'm',
				 'fill_value': i2_nan,
				 'min_value': -10000,
				 'max_value': 10000,
				 'type': 'i2',
				 'scale': 1e-4
				 }
longname = 'Correction processed from KaRIn/NadirAlt cross-overs'
karin_na_height_corr = {'varname': 'KaRIn_NA_height_corr',
                        'longname': longname,
			'unit': 'm',
			'fill_value': i2_nan,
			'min_value': -10000,
			'max_value': 10000,
			'type': 'i2',
			'scale': 1e-4
			}
karin_cal_flag = {'varname': 'KaRIn_cal_flag',
                  'longname': 'quality of the cross-over calibrations',
		  'unit': 'bits',
		  'fill_value': u1_nan,
		  'min_value': 0,
		  'max_value': 255,
		  'type': 'u1',
		  'scale': 1
		  }
