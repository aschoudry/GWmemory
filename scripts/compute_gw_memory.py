import numpy as np
import matplotlib.pyplot as plt
import h5py
from math import factorial as fac

file_name="/home/ashok/gravitational_wave_memory_project/data/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
f = h5py.File(file_name,'r+')  

lmax = 8
time_idx =10


def Y_lm_2_coeff(l, m, time_idx):
	data = f['Extrapolated_N4.dir']['Y_l'+str(l)+'_m'+str(m)+'.dat'][:][time_idx]
	time = data[0]
	Ylm_real_coeff = data[1]
	Ylm_imag_coeff = data[2]
	return time, Ylm_real_coeff, Ylm_imag_coeff 


def Y_lm_m2(l, m, theta, phi):
	k1=max(0, m-2)
	k2=min(l+m, l-2)
	
	dlm2=0
	for k in range(k1,k2+1):
		dlm2+=pow(-1,k)*np.sqrt(fac(l+m)*fac(l-m)*fac(l+2)*fac(l-2))*pow(np.cos(theta/2.0), 2*l+m-2*k-2)*pow(np.sin(theta/2.0), 2*k-m+2)/ \
			(fac(k)*fac(k-m+2)*fac(l+m-k)*fac(l-k-2))

	Ylm_m2 = np.sqrt(2*l+1)*dlm2*np.exp(1j*m*phi)/np.sqrt(4*np.pi)
	return Ylm_m2

def d2h_dt2(time_idx, theta, phi, lmax):
	sum=0
	for l in range(2,lmax+1):
		m_value=range(-l, l+1)
		for m in m_value:
			sum+=Y_lm_2_coeff(l, m, time_idx)[1]*Y_lm_m2(l, m, theta, phi).real + 1j*Y_lm_2_coeff(l, m, time_idx)[2]*Y_lm_m2(l, m, theta, phi).imag
	return sum

def d2h22_by_dt2(time_idx):
	data =  f['Extrapolated_N4.dir']['Y_l2_m2.dat'][:][time_idx]
	time = data[0]
	Ylm_real_coeff = data[1]
	Ylm_imag_coeff = data[2]
	return time, Ylm_real_coeff, Ylm_imag_coeff 


	
time_idxVec=range(0,5000,100)

time_vec=np.array([])
d2h22_values=np.array([])

for i in range(len(time_idxVec)):
	d2h22_values=np.append(d2h22_values,d2h22_by_dt2(time_idxVec[i])[1])
	time_vec=np.append(time_vec, d2h22_by_dt2(time_idxVec[i])[0])

h22_values = np.cumsum(d2h22_values)

plt.plot(time_vec,d2h22_values)
plt.show()
