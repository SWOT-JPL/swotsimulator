import numpy
import matplotlib.pylab as plt
import time
import myspectools
from math import pi
import swotsimulator.const as const
import params as p
import swotsimulator.rw_data as rw_data


dirname=p.outdatadir
rootname=p.file_output
dx=p.delta_al #km
listfile=glob.glob(p.file_output+'_c*'+'_p*.nc')

nr=0
for ifile in listfile:
    print(ifile)
    hh_roll=[]
    hh_phase=[]
    hh_bd=[]
    hh_timing=[]
    data=rw_data.Sat_SWOT(nfile=ifile)
    data.load_swath(wt=[], wt_err=[])
    nal, nac=numpy.shape(data.roll_err)

    if nal*dx>=distance_min:
      f0=numpy.linspace(1/distance_min*dx, 1/2-1/distance_min*dx,num=int(distance_min/dx/2.))
      tap=0.04
      if p.roll: ff, PSD_roll=myspectools.psd1d(hh=data.roll_err,dx=dx, detrend=True, tap=tap)
      if p.phase: ff, PSD_phase=myspectools.psd1d(hh=data.phase_err,dx=dx, detrend=True, tap=tap)
      if p.baselinedilation: ff, PSD_bd=myspectools.psd1d(hh=data.bd_err,dx=dx, detrend=True, tap=tap)
      if p.timing: ff, PSD_timing=myspectools.psd1d(hh=data.timing_err,dx=dx, detrend=True, tap=tap)
      try:

        SS_roll=SS_roll+numpy.interp(f0, ff, PSD_roll)
        SS_phase=SS_phase+numpy.interp(f0, ff, PSD_phase)
        SS_bd=SS_bd+numpy.interp(f0, ff, PSD_bd)
        SS_timing=SS_timing+numpy.interp(f0, ff, PSD_timing)
      except:
        SS_roll=numpy.interp(f0, ff, PSD_roll)
        SS_phase=numpy.interp(f0,ff, PSD_phase)
        SS_bd=numpy.interp(f0, ff, PSD_bd)
        SS_timing=numpy.interp(f0, ff, PSD_timing)
      nr+=1

SS_roll/=nr
SS_phase/=nr
SS_bd/=nr
SS_timing/=nr

## -- COnvert spectrum to be compared with Esteban et al:  
MSH_roll=SS_roll*(1e-6*37.870*1e3)**2
PSHroll=PSroll*(1e-6*37.870*1e3)**2
MSH_phase=SS_phase*(1e-6*37.870*1e3)**2
PSHphase=PSphase*(1e-6*37.870*1e3)**2
MSH_bd=SS_bd*(42.017**2)**2
PSHbd=PSbd*(42.017**2)**2
MSH_timing=SS_timing
PSHtiming=PStiming

MSH_tot=MSH_roll+MSH_phase+MSH_bd+MSH_timing
PSHtot=PSHroll+PSHphase+PSHbd+PSHtiming

ind=numpy.where(((ff>1./40000)&(ff<1./4)))
ff=f0
dff=ff[1]-ff[0]
indi=numpy.where(((spatial_frequency>1./40000)&(spatial_frequency<1./4)))

#IrollPSD1=numpy.sum(rollPSD1[indi]*numpy.diff(spatial_frequency)[indi])
#IrollPSD2=numpy.sum(rollPSD2[indi]*numpy.diff(spatial_frequency)[indi])
#IMS_roll=numpy.sum(MS_roll[ind]/(1+const.sat_elev/const.Rearth)**2*dff)
#IPSHroll=numpy.sum(PSHroll[indi]*1e4*numpy.diff(spatial_frequency)[indi])

plt.close()
plt.loglog(ff,MS_roll/(1+const.sat_elev/const.Rearth)**2, color='k'); plt.grid()
f1=plt.loglog(spatial_frequency,rollPSD1, color='green',lw=2, label='roll knowledge error')
plt.loglog(spatial_frequency,rollPSD2, color='blue',lw=2, label='gyro')
plt.loglog(spatial_frequency,rollPSD1+rollPSD2, color='red',lw=2, label='total roll error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-6,0.25,1e-4,1e7])
plt.xlabel(u'cy/km')
plt.ylabel(u'μrad²/(cy/km)')
plt.legend()
#lt.annotate('integration: '+str(numpy.sqrt(IrollPSD1+IrollPSD2))[:4]+u' μrad' +' --> '+str(numpy.sqrt(IPSHroll))[:4]+u' cm on SSH' , xy=(0.3, 0.75), xycoords='axes fraction', size=10, color='k',weight='bold')
plt.annotate ('', (1./40000,1e4 ), (0.25, 1e4), arrowprops={'arrowstyle':'<->'})
plt.title('Roll error spectra')
plt.savefig('roll_spectra.png')
plt.show()

plt.close()
plt.loglog(ff,MS_phase/((pi/180*1e6)**2 *((1/(const.Fka*2*pi/const.C*const.B)*(1+const.sat_elev/const.Rearth)))**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,phasePSD, color='red',lw=2, label='differential phase drift error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-6,0.25,1e-4,1e8])
plt.xlabel(u'cy/km')
plt.ylabel(u'deg²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,5e5 ), (0.25, 5e5), arrowprops={'arrowstyle':'<->'})
plt.title('Phase error spectra')
plt.savefig('phase_spectra.png')
plt.show()

plt.close()
plt.loglog(ff,MS_bd/(((10**-6)*(1+const.sat_elev/const.Rearth)/(const.sat_elev*const.B)*1e6)**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,dilationPSD, color='red',lw=2, label='baseline dilation error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-6,0.25,1e0,1e7])
plt.xlabel(u'cy/km')
plt.ylabel(u'μm²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,5e5 ), (0.25, 5e5), arrowprops={'arrowstyle':'<->'})
plt.title('Baseline dilation error spectra')
plt.savefig('bd_spectra.png')
plt.show()


plt.close()
plt.loglog(ff,MS_timing/((const.C/2*10**-12)**2), color='k'); plt.grid()
plt.loglog(spatial_frequency,timingPSD, color='red',lw=2, label='timing error')
plt.axvline(1./40000, color='k', linestyle='--')
plt.axis([5e-6,0.25,1e1,1e10])
plt.xlabel(u'cy/km')
plt.ylabel(u'psec²/(cy/km)')
plt.legend()
plt.annotate ('', (1./40000,1e8 ), (0.25, 1e8), arrowprops={'arrowstyle':'<->'})
plt.title('Timing error spectra')
plt.savefig('timing_spectra.png')
plt.show()


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
plt.axis([5e-6,0.25,1e-3,1e8])
plt.xlabel(u'cy/km')
plt.ylabel(u'cm2/(cy/km)')
plt.legend()
plt.title('Total Error spectra')
plt.savefig('{}_allssh_spectra.png'.format(p.config))
plt.show()

