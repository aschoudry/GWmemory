import numpy as np
import matplotlib.pyplot as plt
import plotsettings

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

file_Spinning_binary_with_SpinAlignedN4 = ["rMPsi4_Sz1_0p20_Sz2_0p20_q1p5dataN4", "rMPsi4_Sz1_0p60_Sz2_0p60_q1p5dataN4", "rMPsi4_Sz1_0p80_Sz2_0p80_q1p5dataN4", "rMPsi4_Sz1_0p99_Sz2_0p99_q1p5dataN4"]
file_Spinning_binary_with_SpinAlignedEout = ["rMPsi4_Sz1_0p20_Sz2_0p20_q1p5dataEout", "rMPsi4_Sz1_0p60_Sz2_0p60_q1p5dataEout", "rMPsi4_Sz1_0p80_Sz2_0p80_q1p5dataEout", "rMPsi4_Sz1_0p99_Sz2_0p99_q1p5dataEout"]

file_Spinning_binary_with_SpinAntialignedN4= ["rMPsi4_Sz1_m0p20_Sz2_m0p20_q1p5dataN4", "rMPsi4_Sz1_m0p43_Sz2_m0p43_q1p5dataN4", "rMPsi4_Sz1_m0p60_Sz2_m0p60_q1p5dataN4", "rMPsi4_Sz1_m0p80_Sz2_m0p80_q1p5dataN4", "rMPsi4_Sz1_m0p90_Sz2_m0p90_q1p5dataN4", "rMPsi4_Sz1_m0p94_Sz2_m0p94_q1p5dataN4"]
file_Spinning_binary_with_SpinAntialignedEout= ["rMPsi4_Sz1_m0p20_Sz2_m0p20_q1p5dataEout", "rMPsi4_Sz1_m0p43_Sz2_m0p43_q1p5dataEout", "rMPsi4_Sz1_m0p60_Sz2_m0p60_q1p5dataEout", "rMPsi4_Sz1_m0p80_Sz2_m0p80_q1p5dataEout", "rMPsi4_Sz1_m0p90_Sz2_m0p90_q1p5dataEout", "rMPsi4_Sz1_m0p94_Sz2_m0p94_q1p5dataEout"]


file_Spinning_binary_with_totalSpin0N4 = ["rMPsi4_Sz1_m0p50_Sz2_0p50_q1p5dataN4","rMPsi4_Sz1_m0p60_Sz2_0p60_q1p5dataN4","rMPsi4_Sz1_m0p90_Sz2_0p94_q1p5dataN4"]
file_Spinning_binary_with_totalSpin0Eout = ["rMPsi4_Sz1_m0p50_Sz2_0p50_q1p5dataEout","rMPsi4_Sz1_m0p60_Sz2_0p60_q1p5dataEout","rMPsi4_Sz1_m0p90_Sz2_0p94_q1p5dataEout"]



#import data
#Spinning_binary_with_SpinAligned_27Dec
#Spinning_binary_with_SpinAntialigned_27Dec
#Spinning_binary_with_totalSpin0_27Dec


for filename in file_Spinning_binary_with_totalSpin0Eout:

	#filename = "rMPsi4_noSpin_q1data"
	timeNR, psi_plus, psi_cross = np.loadtxt("/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_totalSpin0_27Dec/"+filename+".txt", unpack=True)

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

	f = open("/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_totalSpin0_27Dec/"+filename+"Clean.txt","w") 

	for i in range(len(timeNR)):
		f.write("%E %E %E\n" % (timeNR[i], psi_plus[i], psi_cross[i])) 
	f.close() 
