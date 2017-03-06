import glob
import numpy
import matplotlib.pylab as plt
import myspectools
from math import pi
import swotsimulator.const as const
import params as p
import swotsimulator.rw_data as rw_data


## -- Load instrumental error -- ## 
file_instr = p.file_inst_error
fid_instr = rw_data.file_instr(file=p.file_inst_error)
fid_instr.read_var(spatial_frequency = [], rollPSD = [], gyroPSD = [],\
                   phasePSD = [], dilationPSD = [], timingPSD = [])
spatial_frequency = fid_instr.spatial_frequency 
rollPSD = fid_instr.rollPSD
gyroPSD = fid_instr.gyroPSD
phasePSD = fid_instr.phasePSD
dilationPSD = fid_instr.dilationPSD
timingPSD = fid_instr.timingPSD

## Roll in microrad**2
PSroll = (rollPSD+gyroPSD)*(1+const.sat_elev/const.Rearth)**2*(1./3600*pi/180*1e6)**2 
## Phase in micrarad**2
PSphase = phasePSD*(pi/180*1e6)**2 *((1/(const.Fka*2*pi/const.C*const.B)*(1+const.sat_elev/const.Rearth)))**2 # microrad**2
## Baseline dilation in m**-2
PSbd = dilationPSD*((10**-6)*(1+const.sat_elev/const.Rearth)/(const.sat_elev*const.B)*1e6)**2  #Timing in psec**2
PStiming = timingPSD*(const.C/2*10**-12)**2


dirname = p.outdatadir
rootname = p.file_output
dx = p.delta_al #km
listfile = glob.glob(p.file_output+'_c*'+'_p*.nc')

nr = 0
distance_min = 400
f0 = numpy.linspace(1/distance_min*dx, 1/2-1/distance_min*dx,num = int(distance_min/dx/2.))
for coordfile in listfile:
#    print(coordfile)
    data = rw_data.Sat_SWOT(file=coordfile)
    data.load_swath(roll_err = [], phase_err = [], bd_err = [], timing_err=[], x_ac = [])

    nal, nac = numpy.shape(data.roll_err)
    if nal*dx>=distance_min:
      tap = 0.04
      if p.roll:     
          roll_err = data.roll_err[:,0]/(data.x_ac[0]*1000)*1.e6
          ffr, PSD_roll = myspectools.psd1d(hh=roll_err,dx=dx, detrend=True, tap=tap)
      if p.phase: 
          phase_err = data.phase_err[:,0]/(data.x_ac[0]*1000)*1.e6
          ffp, PSD_phase = myspectools.psd1d(hh=phase_err,dx=dx, detrend=True, tap=tap)
      if p.baseline_dilation: 
          bd_err = data.bd_err[:,-1]/((data.x_ac[-1])**2)
          ffb, PSD_bd = myspectools.psd1d(bd_err,dx=dx, detrend=True, tap=tap)
      if p.timing: 
          timing_err = data.timing_err[:,-1]
          fft, PSD_timing = myspectools.psd1d(hh=timing_err,dx=dx, detrend=True, tap=tap)

      try: 

        SS_roll = SS_roll+numpy.interp(f0, ffr, PSD_roll)
        SS_phase = SS_phase+numpy.interp(f0, ffp, PSD_phase)
        SS_bd = SS_bd+numpy.interp(f0, ffb, PSD_bd)
        SS_timing = SS_timing+numpy.interp(f0, fft, PSD_timing)
      except:
      
        SS_roll = numpy.interp(f0, ffr, PSD_roll)
        SS_phase = numpy.interp(f0,ffp, PSD_phase)
        SS_bd = numpy.interp(f0, ffb, PSD_bd)
        SS_timing = numpy.interp(f0, fft, PSD_timing)
      nr+=1
  
SS_roll/=nr
SS_phase/=nr
SS_bd/=nr
SS_timing/=nr

## -- Convert spectrum to take into account all values along the swath:   
MSH_roll = SS_roll*(1e-6*37.870*1e3)**2
PSHroll = PSroll*(1e-6*37.870*1e3)**2
MSH_phase = SS_phase*(1e-6*37.870*1e3)**2
PSHphase = PSphase*(1e-6*37.870*1e3)**2
MSH_bd = SS_bd*(42.017**2)**2
PSHbd = PSbd*(42.017**2)**2
MSH_timing = SS_timing
PSHtiming = PStiming

MSH_tot = MSH_roll+MSH_phase+MSH_bd+MSH_timing
PSHtot = PSHroll+PSHphase+PSHbd+PSHtiming

## -- Filter frequencies
ff = f0
ind = numpy.where(((ff>1./40000)&(ff<1./4)))
dff = ff[1]-ff[0]
indi = numpy.where(((spatial_frequency>1./40000)&(spatial_frequency<1./4)))
'''
plt.close()
plt.loglog(ff,MSH_roll/(1+const.sat_elev/const.Rearth)**2, color='k'); plt.grid()
f1=plt.loglog(spatial_frequency,PSroll, color='green',lw=2, label='roll knowledge error')
plt.loglog(spatial_frequency,gyroPSD, color='blue',lw=2, label='gyro')
plt.loglog(spatial_frequency,PSroll, color='red',lw=2, label='total roll error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-3,0.25,1e-4,1e2])
plt.xlabel(u'cy/km')
plt.ylabel(u'μrad²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,1e4 ), (0.25, 1e4), arrowprops={'arrowstyle':'<->'})
plt.title('Roll error spectra')
plt.savefig('roll_spectra.png')
plt.show()

plt.close()
plt.loglog(ff,MSH_phase/((pi/180*1e6)**2 *((1/(const.Fka*2*pi/const.C*const.B)*(1+const.sat_elev/const.Rearth)))**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,phasePSD, color='red',lw=2, label='differential phase drift error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-3,0.25,1e-4,1e2])
plt.xlabel(u'cy/km')
plt.ylabel(u'deg²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,5e5 ), (0.25, 5e5), arrowprops={'arrowstyle':'<->'})
plt.title('Phase error spectra')
plt.savefig('phase_spectra.png')
plt.show()

plt.close()
plt.loglog(ff,MSH_bd/(((10**-6)*(1+const.sat_elev/const.Rearth)/(const.sat_elev*const.B)*1e6)**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,dilationPSD, color='red',lw=2, label='baseline dilation error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-3,0.25,1e0,1e2])
plt.xlabel(u'cy/km')
plt.ylabel(u'μm²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,5e5 ), (0.25, 5e5), arrowprops={'arrowstyle':'<->'})
plt.title('Baseline dilation error spectra')
plt.savefig('bd_spectra.png')
plt.show()


plt.close()
plt.loglog(ff,MSH_timing/((const.C/2*10**-12)**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,timingPSD, color='red',lw=2, label='timing error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-3,0.25,1e1,1e2])
plt.xlabel(u'cy/km')
plt.ylabel(u'psec²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,1e8 ), (0.25, 1e8), arrowprops={'arrowstyle':'<->'})
plt.title('Timing error spectra')
plt.savefig('timing_spectra.png')
plt.show()
'''

plt.close()
plt.loglog(ff,MSH_roll*1e4, color='grey'); plt.grid()
plt.loglog(ff,MSH_phase*1e4, color='grey'); plt.grid()
plt.loglog(ff,MSH_bd*1e4, color='grey'); plt.grid()
plt.loglog(ff,MSH_timing*1e4, color='grey'); plt.grid()
plt.loglog(ff,MSH_tot*1e4, color='k'); plt.grid()
plt.loglog(spatial_frequency,PSHtot*1e4, color='red',lw=2, label='total error')
plt.loglog(spatial_frequency,PSHroll*1e4, color='purple',lw=2, label='roll error')
plt.loglog(spatial_frequency,PSHphase*1e4, color='blue',lw=2, label='phase error')
plt.loglog(spatial_frequency,PSHbd*1e4, color='green',lw=2, label='baseline error')
plt.loglog(spatial_frequency,PSHtiming*1e4, color='brown',lw=2, label='timing error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axvline(7.4e-5, color='k', linestyle='--')
plt.axvline(0.00027, color='k', linestyle='--')
plt.axis([5e-3,0.25,1e-3,1e2])
plt.xlabel(u'cy/km')
plt.ylabel(u'cm2/(cy/km)')
plt.legend()
plt.title('Total Error spectra')
plt.savefig('allssh_spectra.png')
plt.show()

