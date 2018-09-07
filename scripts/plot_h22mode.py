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
	h22_real_SXS=np.append(h22_real_SXS,data_h[:][i][1])
	h22_imag_SXS=np.append(h22_imag_SXS,data_h[:][i][2])


#intergrating the Phi4 data with respect to time first time
dhSXS_dt_initial_real=np.gradient(h22_real_SXS, time)[0]
dhSXS_dt_initial_imag=np.gradient(h22_imag_SXS, time)[0]

dh22_real_by_dt=integrate.cumtrapz(d2h22_real_by_dt2, time,initial=0*dhSXS_dt_initial_real)
dh22_imag_by_dt=integrate.cumtrapz(d2h22_imag_by_dt2, time, initial=0*dhSXS_dt_initial_real)

#integrating the data second time with respect to time

h22_real=integrate.cumtrapz(dh22_real_by_dt, time, initial=0*h22_real_SXS[0])
h22_imag=integrate.cumtrapz(dh22_imag_by_dt, time, initial=0*h22_imag_SXS[0])

#making plots
plt.plot(time, h22_real, 'g', label="h22 real")
plt.plot(time, h22_real_SXS,'r--', label="h22 real SXS")
plt.xlabel('time')
plt.ylabel('h22 amplitude')
plt.legend()
plt.show()

plt.plot(time, h22_real-h22_real_SXS)
plt.plot(time, h22_imag-h22_imag_SXS)
plt.show()

#plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/h22_plot2.pdf")

# Plot Psi4 data
dh22_real_SXS_by_dt = np.gradient(h22_real_SXS, time)
dh22_imag_SXS_by_dt = np.gradient(h22_imag_SXS, time)

d2h22_real_SXS_by_dt2=np.gradient(dh22_real_SXS_by_dt, time)
d2h22_imag_SXS_by_dt2=np.gradient(dh22_imag_SXS_by_dt, time)

plt.plot(time,d2h22_real_SXS_by_dt2,'r--' )
plt.plot(time,d2h22_real_by_dt2,'k' )
plt.show()

#Subtract linearfit
linear_fit_coeff_real = np.polyfit(time, h22_real, 1)
linear_fit_coeff_imag = np.polyfit(time, h22_imag, 1)

linear_fit_real=linear_fit_coeff_real[0]*time + linear_fit_coeff_real[1]
plt.plot(time, h22_real-linear_fit_real,'r--', label="linear fit sub from h22")
plt.plot(time, h22_real_SXS, 'g', label="h22 SXS")
plt.legend()
plt.show()  






  
 





