import numpy as np
import matplotlib.pyplot as plt
import plotsettings

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

filename_list_nonspinning = ["rMPsi4_noSpin_q1data"]#, "rMPsi4_noSpin_q1p5data", "rMPsi4_noSpin_q2data", "rMPsi4_noSpin_q2p5data", "rMPsi4_noSpin_q3data", "rMPsi4_noSpin_q4data", "rMPsi4_noSpin_q5data", "rMPsi4_noSpin_q6data", "rMPsi4_noSpin_q7data", "rMPsi4_noSpin_q8data", "rMPsi4_noSpin_q9p5data"]

#filename_list_equalmass_alignSzSpin = ["rMPsi4_AlignedSpin_Sz_0p2_q1data", "rMPsi4_AlignedSpin_Sz_0p43_q1data", "rMPsi4_AlignedSpin_Sz_0p600_q1data", "rMPsi4_AlignedSpin_Sz_0p6_q1data", "rMPsi4_AlignedSpin_Sz_0p800_q1data", "rMPsi4_AlignedSpin_Sz_0p95_q1data", "rMPsi4_AlignedSpin_Sz_m0p2_q1data", "rMPsi4_AlignedSpin_Sz_m0p9_q1data", "rMPsi4_AlignedSpin_Sz_mp800_q1data"]

#import data

for filename in filename_list_nonspinning:

	#filename = "rMPsi4_noSpin_q1data"
	timeNR, psi_plus, psi_cross = np.loadtxt("../data/NonSpinning_differentMassRatio/"+filename+".txt", unpack=True)

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

	f = open("../data/NonSpinning_differentMassRatio/"+filename+"Clean.txt","w") 

	for i in range(len(timeNR)):
		f.write("%E %E %E\n" % (timeNR[i], psi_plus[i], psi_cross[i])) 
	f.close() 
