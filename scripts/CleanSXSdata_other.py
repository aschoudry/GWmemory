import numpy as np
import matplotlib.pyplot as plt
import plotsettings

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

filename_list_spinning_algned_with_different_mag = ["rMPsi4_Sz1_m0p4_Sz2_0p8_q1p5data", "rMPsi4_Sz1_m0p5_Sz2_m0p9_q1p5data", "rMPsi4_Sz1_m0p625_Sz2_m0p25_q1p5data", "rMPsi4_Sz1_m0p9_Sz2_m0p5_q1p5data"]

filename_list_spinning_only_one_Sz = ["rMPsi4_Sz1_m0p3_q1p5data", "rMPsi4_Sz1_m0p6_q1p5data", "rMPsi4_Sz1_m0p9_q1p5data"]
filename_list_spinning_Sz_antialigned = ["rMPsi4_Sz1_m0p5_Sz2_0p5_q1p5data", "rMPsi4_Sz1_m0p9_Sz2_0p9_q1p5data"]

#import data

for filename in filename_list_spinning_Sz_antialigned:

	#filename = "rMPsi4_noSpin_q1data"
	timeNR, psi_plus, psi_cross = np.loadtxt("/home/ashok/gravitational_wave_memory_project/data/Spinning_Sz_antialigned/"+filename+".txt", unpack=True)

	plt.plot(timeNR, psi_plus)
	plt.plot(timeNR, psi_cross)
	plt.show()

	psi_tot=psi_plus + 1j*psi_cross
	phase=np.unwrap(np.angle(psi_tot))
	instfreq=np.gradient(phase, timeNR)

	plt.plot(timeNR, instfreq)
	plt.show()

	lower_cutoff_time=input("enter cutoff")
	higher_cutoff_time=input("hiher cutoff")

	idx_ini =find_nearest_idx(timeNR, lower_cutoff_time)
	idx_fin =find_nearest_idx(timeNR, higher_cutoff_time)

	timeNR=timeNR[idx_ini:idx_fin]
	instfreq=instfreq[idx_ini:idx_fin]
	psi_plus=psi_plus[idx_ini:idx_fin]
	psi_cross=psi_cross[idx_ini:idx_fin]

	plt.plot(timeNR,instfreq)
	plt.show()

	f = open("/home/ashok/gravitational_wave_memory_project/data/Spinning_Sz_antialigned/"+filename+"Clean.txt","w") 

	for i in range(len(timeNR)):
		f.write("%E %E %E\n" % (timeNR[i], psi_plus[i], psi_cross[i])) 
	f.close() 
