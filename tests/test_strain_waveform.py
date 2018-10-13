import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import h5py

# Reading the data file for Phi4
file_name_psi="/home/ashok/gravitational_wave_memory_project/data/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
f_psi = h5py.File(file_name_psi,'r+')

# Reading the strain data
file_name_h="/home/ashok/gravitational_wave_memory_project/data/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
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
	h22_real_SXS=np.append(h22_real_SXS,data_h[:][i][2])
	h22_imag_SXS=np.append(h22_imag_SXS,data_h[:][i][2])

plt.plot(time, h22_real_SXS)
plt.show()
