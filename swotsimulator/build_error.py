# import params as p
import datetime
import numpy
import swotsimulator.rw_data as rw_data
import swotsimulator.const as const
import swotsimulator.mod_tools as mod_tools
import swotsimulator.error as comp_error
from math import pi, sqrt
import netCDF4
import sys
from scipy.ndimage.filters import gaussian_filter
import logging

# Define logging level for debug purposes
logger = logging.getLogger(__name__)


def reconstruct_2D_error(x_ac, err_out, dict_noise):
    nac = numpy.shape(x_ac)[0]
    ncenter = int(nac / 2)
    ac_l = numpy.mat(x_ac[:ncenter])
    ac_r = numpy.mat(x_ac[ncenter:])
    if 'phase' in dict_noise.keys():
        phase1d = dict_noise['phase']
        err_out.phase[:, :ncenter] = numpy.mat(phase1d[:, 0]).T * ac_l
        err_out.phase[:, ncenter:] = numpy.mat(phase1d[:, 1]).T * ac_r
    if 'roll' in dict_noise.keys():
        ac = numpy.mat(x_ac)
        roll1d = dict_noise['roll']
        err_out.roll[:, :] = numpy.mat(roll1d).T * ac
    if 'corrected_roll_phase' in dict_noise.keys():
        rollphase_est_1d = dict_noise['corrected_roll_phase']
        rollphase_est = numpy.full((rollphase_est_1d.shape[0], nac), numpy.nan)
        rollphase_est[:, :ncenter] = numpy.mat(rollphase_est_1d[:,0]).T * ac_l
        rollphase_est[:, ncenter:] = numpy.mat(rollphase_est_1d[:, 1]).T * ac_r
        err_out.corrected_roll_phase[:, :] = rollphase_est
    if 'baseline_dilation' in dict_noise.keys():
        ac2 = numpy.mat((x_ac)**2)
        baseline_dilation1d = numpy.mat(dict_noise['baseline_dilation'])
        err_out.baseline_dilation[:, :] = baseline_dilation1d.T * ac2
    if 'timing' in dict_noise.keys():
        timing1d = dict_noise['timing']
        ones_ac = numpy.mat(numpy.ones((int(nac/2))))
        err_out.timing[:, :int(nac / 2)] = numpy.mat(timing1d[:, 0]).T*ones_ac
        err_out.timing[:, int(nac / 2):] = numpy.mat(timing1d[:, 1]).T*ones_ac
    return None



class error():
    '''Class error define all the possible errors that can be computed using
    SWOT simulator.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in
    file file_coeff, the random realisations are read directly using
    load_coeff.  The corresponding errors on a swath can be computed using
    make_error. '''
    def __init__(self, p, roll=None, ssb=None, wet_tropo=None, phase=None,
                 baseline_dilation=None, karin=None, timing=None, SSH=None,
                 wt=None):
        self.roll = roll
        self.ssb = ssb
        self.wet_tropo = wet_tropo
        self.phase = phase
        self.baseline_dilation = baseline_dilation
        self.karin = karin
        self.timing = timing
        self.SSH_error = SSH
        self.wt = wt
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand

    def init_error(self, p, nac, seed=0):
        '''Initialization of errors: Random realisation of errors are
        computed using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of
        each random realisation.
        By default, there are ncomp1d=2000 random realisations for
        the instrumental errors (1d spectrum) and ncomp2d=2000 random
        realisations for the geophysical errors (2d spectrum)
        and nrandkarin*x_ac km of random number for KaRIN noise.'''
        numpy.random.seed(seed)
        if p.error_spectrum is not None:
            if p.lambda_cut is None:
                p.lambda_cut = 20000
            file_instr = rw_data.file_instr(file=p.error_spectrum)
            spatial_frequency = []
            file_instr.read_var(spatial_frequency=spatial_frequency)
            # - Cut frequencies larger than Nyquist frequency and cut long
            #   wavelengths (larger than p.lambda_max)
            cut_min = 1. / float(2 * p.delta_al)
            cut_max = 1. / p.lambda_max
            ind = numpy.where((file_instr.spatial_frequency < cut_min)
                              & (file_instr.spatial_frequency > 0)
                              & (file_instr.spatial_frequency > cut_max))[0]
            freq = file_instr.spatial_frequency[ind]
            f_cut = 1. / p.lambda_cut
            ind_cut = numpy.where(freq < f_cut)
            min_ind_cut = numpy.min(numpy.where(freq >= f_cut))
            self.spatial_frequency = freq
            if 'RollPhase' in p.noise:
                # - Read and compute roll power spectrum, wavelength longer
                #   than p.lambda_cut are supposed to be corrected
                PSroll = []
                PSgyro = []
                file_instr.read_var(rollPSD=PSroll, gyroPSD=PSgyro)
                PSTroll = file_instr.rollPSD[ind]
                PSgyro = file_instr.gyroPSD[ind]
                PSTroll[ind_cut] = PSTroll[min_ind_cut]
                PSgyro[ind_cut] = PSgyro[min_ind_cut]
                self.roll_psd = PSTroll
                self.gyro_psd = PSgyro
            if 'RollPhase' in p.noise:
                # - Read Phase power spectrum, wavelength longer than
                #   p.lambda_cut are supposed to be corrected
                phasePSD = []
                file_instr.read_var(phasePSD=phasePSD)
                PSphase = file_instr.phasePSD[ind]
                PSphase[ind_cut] = PSphase[min_ind_cut]
                self.phase_psd = PSphase
                # - Compute left and right random coefficients using
                #   the power spectrum
            if 'BaselineDilation' in p.noise:
                # - Read baseline dilation power spectrum, wavelength longer
                #   than p.lambda_cut are supposed to be corrected
                dilationPSD = []
                file_instr.read_var(dilationPSD=dilationPSD)
                PSbd = file_instr.dilationPSD[ind]
                PSbd[ind_cut] = PSbd[min_ind_cut]
                self.dilation_psd = PSbd
            if 'Timing' in p.noise:
                # - Read timing power spectrum, wavelength longer than
                #   p.lambda_cut are supposed to be corrected
                timingPSD = []
                file_instr.read_var(timingPSD=timingPSD)
                PStim = file_instr.timingPSD[ind]
                PStim[ind_cut] = PStim[min_ind_cut]
                self.timing_psd = PStim
            if 'WetTroposphere' in p.noise:
                # - Define power spectrum of error in path delay
                #   due to wet tropo
                f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
                # - Global mean wet tropo power spectrum in
                #   cm**2/(cycle/km) for L >= 100 km
                PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
                # - Wet tropo power spectrum in cm**2/(cycle/km)
                #   for L < 100 km
                indf = numpy.where(f > 10**-2)
                PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
                # - Compute random coefficients in 2D using the previously
                #   defined power spectrum
                gencoef = mod_tools.gen_coeff_signal2d(f, PSwt, self.ncomp2d)
                self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
                # - Define radiometer error power spectrum for a beam
                #   High frequencies are cut to filter the associated error:
                #   during the reconstruction of the wet trop signal
                # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
                PSradio = 9.5 * 10**-5 * f**-1.79
                PSradio[numpy.where((f < 1./1000.))] = \
                    9.5 * 10**-5 * (10**-3)**-1.79
                indf = numpy.where((f > 0.0023) & (f <= 0.0683))
                PSradio[indf] = 0.036 * f[indf]**-0.814
                PSradio[numpy.where(f > 0.0683)] = 0.32
                # - Compute random coefficients (1D) for the radiometer error
                #   power spectrum for right and left beams
                gencoef = mod_tools.gen_coeff_signal1d(f, PSradio,
                                                       self.ncomp2d)
                self.A_radio_r, self.phi_radio_r, self.fr_radio_r = gencoef
                self.A_radio_l, self.phi_radio_l, self.fr_radio_l = gencoef


    def make_error(self, sgrid, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, the phase between the two signals,
        the timing error, the roll of the satellite, the sea surface bias,
        the distorsion of the mast,
        the karin noise due to the sensor itself. '''
        # - Errors array are the same size as the swath size
        nal, nac = numpy.shape(SSH_true)
        x_al = sgrid.x_al + sgrid.al_cycle * cycle
        # ind_al=numpy.arange(0,nal)
        if 'Karin' in p.noise:
            error_karin = comp_error.Karin(p)
            # TODO tmp swh varying in space
            swh = 0 * SSH_true + p.swh
            seed = int(sgrid.x_al[0]+sgrid.al_cycle)
            dic_error = error_karin.generate(seed, sgrid.x_al, sgrid.x_ac,
                                              swh)
            self.karin = dic_error["simulated_error_karin"]
        if 'CorrectedRollPhase' in p.noise:
            first_time = numpy.datetime64(p.first_time)
            rollphase = comp_error.CorrectedRollPhase(p, first_time)
            results = rollphase._generate_1d(sgrid.time)
            self.roll1d, self.phase1d, self.rollphase_est1d = results
        if 'RollPhase' in p.noise:
            rollphase = comp_error.RollPhase(p, self.roll_psd, self.gyro_psd,
                                             self.phase_psd,
                                             self.spatial_frequency)
            self.roll1d, self.phase1d = rollphase._generate_1d(x_al)
        if 'BaselineDilation' in p.noise:

            bd = comp_error.BaselineDilation(p, self.dilation_psd,
                                             self.spatial_frequency)
            self.baseline_dilation1d = bd._generate_1d(x_al)
        if 'Timing' in p.noise:
            timing = comp_error.Timing(p, self.timing_psd,
                                       self.spatial_frequency)
            self.timing1d = timing._generate_1d(x_al)
        if 'WetTroposphere' in p.noise:
            wt = comp_error.WetTroposphere(p)
            dic_error = wt.generate(x_al, sgrid.x_ac)
            self.wet_tropo2 = dic_error["simulated_error_troposphere"]
            self.wet_tropo2nadir = dic_error["simulated_error_troposphere_nadir"]
            # self.wtnadir
        #if p.ssb is True:
        #    logger.error("No SSB error implemented yet")
        return None

    def reconstruct_2D(self, p, x_ac):
        dict_noise = {}
        if 'RollPhase' in p.noise:
            dict_noise['phase'] = self.phase1d
            dict_noise['roll'] = self.roll1d
        if 'CorrectedRollPhase' in p.noise:
            dict_noise['phase'] = self.phase1d
            dict_noise['roll'] = self.roll1d
            dict_noise['corrected_roll_phase'] = self.rollphase_est1d
        if 'BaselineDilation' in p.noise:
            dict_noise['baseline_dilation'] = self.baseline_dilation1d
        if 'Timing' in p.noise:
            dict_noise['timing'] = self.timing1d
        reconstruct_2D_error(x_ac, self, dict_noise)
        return None

    def make_SSH_error(self, SSH_true, p):
        '''Compute observed SSH adding all the computed error to the model SSH.
        If residual path delay errors after 2-beams and 1-beam radiometer
        correction are both computed (nbeam='both'), only the path delay error
        after 1-beam radiometer correction is considered in the
        observed SSH. '''
        numpy.seterr(invalid='ignore')
        self.SSH = SSH_true
        if 'Karin' in p.noise:
            self.SSH = self.SSH + self.karin
        if 'Timing' in p.noise:
            self.SSH = self.SSH + self.timing
        if 'CorrectedRollPhase' in p.noise:
            self.SSH = self.SSH + self.corrected_roll_phase
        elif 'RollPhase' in p.noise:
            self.SSH = self.SSH + self.roll
            self.SSH = self.SSH + self.phase
        if 'BaselineDilation' in p.noise:
            self.SSH = self.SSH + self.baseline_dilation
        if 'WetTroposphere' in p.noise:
                self.SSH = self.SSH + self.wet_tropo2
        if p.file_input is not None:
            self.SSH[numpy.where(SSH_true == p.model_nan)] = p.model_nan



class errornadir():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in file
    file_coeff, the random realisations are read directly using load_coeff.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self, p, nadir=None, wet_tropo1=None, wt=None,
                 nadir_alone=False):
        self.nadir = nadir
        self.wet_tropo1 = wet_tropo1
        self.wt = wt
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d
        self.nbeam = 1
        self.nadir_alone = nadir_alone

    def init_error(self, p, seed=0, nadir_alone=False):
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        numpy.random.seed(seed)
        # Run reprodictible: generate or load nrand random numbers:
        # - Compute random coefficients in 1D for the nadir error
        wnoise = getattr(p, 'wnoise', 100)
        p.wnoise = wnoise
        # - Define the sepctrum of the nadir instrument error
        # self.A=numpy.random.normal(0.0,sqrt(p.wnoise)
        # /numpy.float64(sqrt(2*p.delta_al)), (self.nrand))*0.01
        f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = numpy.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        PSD = PSD * 10**(-4)
        self.freq = f
        self.psd = PSD
        if 'WetTroposphere' in p.noise and self.nadir_alone == True:
            # - Define power spectrum of error in path delay due to wet tropo
            f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
            # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
            #   for L >= 100 km
            PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
            # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
            indf = numpy.where(f > 10**-2)
            PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
            # - Compute random coefficients in 2D using the previously defined
            #   power spectrum
            gencoef = mod_tools.gen_coeff_signal2d(f, PSwt, self.ncomp2d)
            self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
            # - Define radiometer error power spectrum for a beam
            #   High frequencies are cut to filter the associated error during
            #   the reconstruction of the wet trop signal
            # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
            PSradio = 9.5 * 10**-5 * f**-1.79
            PSradio[numpy.where((f < 1./1000.))] = 9.5 * 10**-5*(10**-3)**-1.79
            indf = numpy.where((f > 0.0023) & (f <= 0.0683))
            PSradio[indf] = 0.036 * f[indf]**-0.814
            PSradio[numpy.where(f > 0.0683)] = 0.32
            # - Compute random coefficients (1D) for the radiometer error power
            #   spectrum for right and left beams
            gencoef = mod_tools.gen_coeff_signal1d(f, PSradio, self.ncomp2d)
            self.A_radio, self.phi_radio, self.fr_radio = gencoef
        return None


    def make_error(self, orb, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        nal = numpy.shape(SSH_true)[0]
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        # - Compute the correspond timing error on the swath in m
        x_al = orb.x_al + orb.al_cycle * cycle
        if 'Altimeter' in p.noise:
            alt = comp_error.Altimeter(p)
            dic_error = alt.generate(x_al)
            self.nadir = dic_error["simulated_error_altimeter"]
        if 'WetTroposphere' in p.noise and self.nadir_alone == True:
            wt = comp_error.WetTroposphere(p)
            dic_error = wt.generate(x_al, numpy.arange(-1, 2, 1))
            self.wet_tropo1 = dic_error["simulated_error_troposphere"]
            self.wet_tropo1 = dic_error["simulated_error_troposphere_nadir"]
            self.wet_tropo1nadir = dic_error["simulated_error_troposphere_nadir"]
        return None
