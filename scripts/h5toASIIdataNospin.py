import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import h5py

filename_list_noSpin = ["rMPsi4_noSpin_q1data"]#, "rMPsi4_noSpin_q1p5data", "rMPsi4_noSpin_q2p5data", "rMPsi4_noSpin_q4data","rMPsi4_noSpin_q6data", "rMPsi4_noSpin_q8data", "rMPsi4_noSpin_q2data", "rMPsi4_noSpin_q3data", "rMPsi4_noSpin_q5data", "rMPsi4_noSpin_q7data", "rMPsi4_noSpin_q9p5data"]

#filename_list_AlignedSpin = ["rMPsi4_AlignedSpin_Sz_0p2_q1data", "rMPsi4_AlignedSpin_Sz_0p43_q1data", "rMPsi4_AlignedSpin_Sz_0p600_q1data", "rMPsi4_AlignedSpin_Sz_0p6_q1data", "rMPsi4_AlignedSpin_Sz_0p800_q1data", "rMPsi4_AlignedSpin_Sz_0p95_q1data", "rMPsi4_AlignedSpin_Sz_m0p2_q1data", "rMPsi4_AlignedSpin_Sz_m0p9_q1data", "rMPsi4_AlignedSpin_Sz_mp800_q1data"]

#filename = "rMPsi4_noSpin_q1data"
# Reading the data file for Phi4

for filename in filename_list_noSpin: 

	file_name_psi="/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/"+filename+".h5"
	f_psi = h5py.File(file_name_psi,'r+')

	# Reading the strain data
	#file_name_h="/users/aschoudhary/gravitational_wave_memory_project/data/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
	#f_h = h5py.File(file_name_h,'r+')

	#reading data for time and l2m2 mode
	data_psi = f_psi['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
	#data_h = f_h['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
	time = np.array([])
	d2h22_real_by_dt2 = np.array([])
	d2h22_imag_by_dt2 = np.array([])
	#h22_real_SXS=np.array([])
	#h22_imag_SXS=np.array([])


	for i in range(len(data_psi)):
		time=np.append(time, data_psi[:][i][0])
		d2h22_real_by_dt2=np.append(d2h22_real_by_dt2, data_psi[:][i][1])
		d2h22_imag_by_dt2=np.append(d2h22_imag_by_dt2, data_psi[:][i][2])
		#h22_real_SXS=np.append(h22_real_SXS,data_h[:][i][1])
		#h22_imag_SXS=np.append(h22_imag_SXS,data_h[:][i][2])

	h22_tot = d2h22_real_by_dt2 + 1j* d2h22_imag_by_dt2
	time_intrp=np.linspace(time[0], time[-1], len(time))
	d2h22_phase=-np.unwrap(np.angle(h22_tot))
	d2h22_amp = abs(h22_tot)

	d2h22_phase_intrp=interp1d(time, d2h22_phase, kind='cubic')
	d2h22_amp_intrp=interp1d(time, d2h22_amp, kind='cubic')

	d2h22_phase_intrp=d2h22_phase_intrp(time_intrp)
	d2h22_amp_intrp=d2h22_amp_intrp(time_intrp)

	#plt.plot(time_intrp, d2h22_amp_intrp)
	#plt.show()

	d2h22_real_by_dt2_intrp_data= (d2h22_amp_intrp*np.exp(1j*d2h22_phase_intrp)).real
	d2h22_imag_by_dt2_intrp_data= (d2h22_amp_intrp*np.exp(1j*d2h22_phase_intrp)).imag


	f = open("/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/"+filename+".txt","w") 

	for i in range(len(time)):
		f.write("%E %E %E\n" % (time_intrp[i], d2h22_real_by_dt2_intrp_data[i], d2h22_imag_by_dt2_intrp_data[i])) 
	f.close() 
	print (filename, "written")
