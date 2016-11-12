''' Roll and gyro power spectrum as a function of wavenumber''' 
import matplotlib.pylab as plt
import numpy
import matplotlib as mpl
import swotsimulator.rw_data as rw_data
import params as p

file_instr = rw_data.file_instr(file=p.file_inst_error)
file_instr.read_var(spatial_frequency=[], phasePSD=[])
## - Cut frequencies larger than Nyquist   frequency and cut long wavelengths (larger than p.lambda_max)
ind=numpy.where((file_instr.spatial_frequency<1./float(2*p.delta_al)) & (file_instr.spatial_frequency>0) & (file_instr.spatial_frequency>1./p.lambda_max))[0]
freq=file_instr.spatial_frequency[ind]
freq2=file_instr.spatial_frequency
fig = plt.figure(figsize=(12,9))
tloc=0.11
tfont=20
stitle = 'Phase power spectrum'
print(freq, file_instr.rollPSD)
#plt.loglog(freq, file_instr.rollPSD[ind], label='roll')
#plt.loglog(freq, file_instr.gyroPSD[ind], label='Gyro')
plt.loglog(freq2, file_instr.phasePSD)
plt.ylabel('Power(asec**2/(cy/km))')
plt.xlabel('Wavenumber (cy/km)')
plt.legend()
plt.grid()
plt.title(stitle, y=-tloc, fontsize=tfont) #size[1])
plt.savefig('Fig9.png')
