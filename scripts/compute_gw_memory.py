import numpy as np
from scipy import integrate
from scipy import signal
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import h5py
from numpy import fft

def compute_integration_usingFFI(dt, w_cut, f):
	
	
	#Compute Fourier transform by numpy's FFT function
	g=np.fft.fft(f)
	#frequency normalization factor is 2*np.pi/dt
	w = np.fft.fftfreq(f.size)*2*np.pi/dt
	
	#In order to get a discretisation of the continuous Fourier transform
	#we need to multiply g by a phase factor
	g*=dt/(np.sqrt(2*np.pi))

	w[0]=pow(10,-30)
	dw=w[2]-w[1]

	for i in range(len(w)):
		if abs(w[i])<=w_cut:
			w[i]=np.sign(w[i])*w_cut

	invFT=np.fft.ifft(-1j*g/w)*len(g)*dw/np.sqrt(2.0*np.pi)
	
	return invFT

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def compute_memory(h_vec, dt):
	mem=np.cumsum(abs(h_vec*h_vec))*dt
	return mem
	

# Reading the data file for Phi4

file_name_psi="/home/ashok/Downloads/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
#file_name_psi="/home/ashok/gravitational_wave_memory_project/data/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
f_psi = h5py.File(file_name_psi,'r+')

# Reading the strain data
file_name_h="/home/ashok/Downloads/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
#file_name_h="/home/ashok/gravitational_wave_memory_project/data/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
f_h = h5py.File(file_name_h,'r+')

#reading data for time and l2m2 mode
data_psi = f_psi['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
data_h = f_h['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
time = np.array([])
d2h22_real_by_dt2 = np.array([])
d2h22_imag_by_dt2 = np.array([])
h22_real_SXS=np.array([])
h22_imag_SXS=np.array([])

for i in range(len(data_psi)):
	time=np.append(time, data_psi[:][i][0])
	d2h22_real_by_dt2=np.append(d2h22_real_by_dt2, data_psi[:][i][1])
	d2h22_imag_by_dt2=np.append(d2h22_imag_by_dt2, data_psi[:][i][2])
	h22_real_SXS=np.append(h22_real_SXS,data_h[:][i][1])
	h22_imag_SXS=np.append(h22_imag_SXS,data_h[:][i][2])

plt.plot(time, d2h22_real_by_dt2)
plt.plot(time, d2h22_imag_by_dt2)
plt.show()

psi_tot=d2h22_real_by_dt2 + 1j*d2h22_imag_by_dt2
phase=-np.unwrap(np.angle(psi_tot))
instfreq=np.gradient(phase, time)

plt.plot(time, instfreq)
plt.show()

plt.plot(time)
plt.show()
idx_ini =find_nearest_idx(time, 100.0)
idx_fin =find_nearest_idx(time, 6230.0)

print idx_ini, idx_fin
time = time[idx_ini:idx_fin]
d2h22_real_by_dt2=d2h22_real_by_dt2[idx_ini:idx_fin]
d2h22_imag_by_dt2=d2h22_imag_by_dt2[idx_ini:idx_fin]
h22_real_SXS=h22_real_SXS[idx_ini:idx_fin]
h22_imag_SXS=h22_imag_SXS[idx_ini:idx_fin]



h22_tot = d2h22_real_by_dt2 + 1j* d2h22_imag_by_dt2
time_intrp=np.linspace(time[0], time[-1], len(time))
d2h22_phase=-np.unwrap(np.angle(h22_tot))
d2h22_amp = abs(h22_tot)

d2h22_phase_intrp=interp1d(time, d2h22_phase)
d2h22_amp_intrp=interp1d(time, d2h22_amp)

d2h22_phase_intrp=d2h22_phase_intrp(time_intrp)
d2h22_amp_intrp=d2h22_amp_intrp(time_intrp)

plt.plot(time_intrp, d2h22_amp_intrp)
plt.show()

d2h22_real_by_dt2_intrp= (d2h22_amp_intrp*np.exp(1j*d2h22_phase_intrp)).real
d2h22_imag_by_dt2_intrp= (d2h22_amp_intrp*np.exp(1j*d2h22_phase_intrp)).imag



plt.plot(time_intrp, d2h22_real_by_dt2_intrp)
plt.show()



#compute_integration_usingFFI(t0, tf, dt, w_cut, f)
t0=time_intrp[0]
tf=time_intrp[-1]
dt=time_intrp[6]-time_intrp[5]
w_cut=0.031
ddf=d2h22_real_by_dt2_intrp


df = compute_integration_usingFFI(dt, w_cut, ddf)
f = compute_integration_usingFFI(dt, w_cut, df)


plt.plot(time_intrp, f, 'k')
plt.plot(time,h22_real_SXS, 'r--' )
plt.show()


mem_hp = compute_memory(df, dt)

plt.plot(time_intrp, mem_hp)
plt.show()
  
 





