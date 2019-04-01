import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
#file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'

#plots setting 
legend_size = 2
fig = plt.figure()
plt.title('Memory Calculation Vs frequency for antialigned spins' ,fontsize = 15)
fontP = FontProperties()
fontP.set_size('20.')
legend = plt.legend(loc='best',prop={'size':legend_size})

#data
#datafile_hMemNR="rMPsi4_Sz1_m0p94_Sz2_m0p94_q1p5dataN4Clean_hMemNR.dat"
#datafile_hphcNR="rMPsi4_Sz1_m0p94_Sz2_m0p94_q1p5dataN4Clean_hNR.dat"

Spin_vec = [-0.2, -0.43, -0.94]
filename_vec=['m0p20', 'm0p43', 'm0p94']

i=0

for filename in filename_vec:

	datafile_hMemNR="rMPsi4_Sz1_"+filename+"_Sz2_"+filename+"_q1p5dataN4Clean_hMemNR.dat"
	datafile_hphcNR="rMPsi4_Sz1_"+filename+"_Sz2_"+filename+"_q1p5dataN4Clean_hNR.dat"


	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	timeNR, hp, hc = np.loadtxt(file_location+datafile_hphcNR, unpack=True)

	h_tot = hp+1j*hc

	phase = np.unwrap(np.angle(h_tot))

	dt = timeNR[2]-timeNR[1]
	phase_dot = np.diff(phase)/dt
	freq = phase_dot/np.pi

	hmem=hmem[1:]
	x=(phase_dot)**(2.0/3.0)

	print len(x), len(hmem)

	plt.plot(x, hmem, label=Spin_vec[i])
	i+=1
#plt.show()

	#plt.plot(timeNR[1:], hmem)

file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'

Spin_vec = [0.2, 0.60, 0.80, 0.99]
filename_vec=['0p20', '0p60', '0p80', '0p99']

i=0

for filename in filename_vec:

	datafile_hMemNR="rMPsi4_Sz1_"+filename+"_Sz2_"+filename+"_q1p5dataN4Clean_hMemNR.dat"
	datafile_hphcNR="rMPsi4_Sz1_"+filename+"_Sz2_"+filename+"_q1p5dataN4Clean_hNR.dat"


	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	timeNR, hp, hc = np.loadtxt(file_location+datafile_hphcNR, unpack=True)

	h_tot = hp+1j*hc

	phase = np.unwrap(np.angle(h_tot))

	dt = timeNR[2]-timeNR[1]
	phase_dot = np.diff(phase)/dt
	freq = phase_dot/np.pi

	hmem=hmem[1:]
	x=(phase_dot)**(2.0/3.0)

	print len(x), len(hmem)

	plt.plot(x, hmem)
	plt.xlim(0.1,0.5)
	plt.ylim(0,0.0009)
	plt.plot(x, hmem, label=Spin_vec[i])
	#plt.show()

	#plt.plot(timeNR[1:], hmem)
	i+=1

plt.legend(loc=2)
plt.savefig("test.pdf")
plt.show()


