import numpy as np
import h5py
from math import factorial as fac

file_name="/home/ashok/gravitational_wave_memory_project/data/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
f = h5py.File(file_name,'r+')  

time_idx =0
data = f['Extrapolated_N2.dir']['Y_l2_m-1.dat'][:][time_idx]


time = data[0]
Ylm_real_coeff = data[1]
Ylm_imag_coeff = data[2]


def Y_lm_m2(l, m, theta, phi):
	k1=max(0, m-2)
	k2=min(l+m, l-2)
	
	dlm2=0
	for k in range(k1,k2+1):
		dlm2+=pow(-1,k)*np.sqrt(fac(l+m)*fac(l-m)*fac(l+2)*fac(l-2))*pow(np.cos(theta/2.0), 2*l+m-2*k-2)*pow(np.sin(theta/2.0), 2*k-m+2)/ \
			(fac(k)*fac(k-m+2)*fac(l+m-k)*fac(l-k-2))

	Ylm_m2 = np.sqrt(2*l+1)*dlm2*np.exp(1j*m*phi)
	return Ylm_m2

def h22(t, theta, phi):
	

