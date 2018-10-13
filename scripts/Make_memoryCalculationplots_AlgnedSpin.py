import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SamMassRatio_differentSzSpinOnly/Memory_data/'
#import data
mass_ratio_vec = [0.43, 0.6, 0.8, 0.95]
filename_vec=['0p43', '0p600', '0p800','0p95']
i=0

for filename in filename_vec:
	
	datafile_hNRdot='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hdotNR.dat'
	datafile_hNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hNR.dat'
	datafile_hMemNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hMemNR.dat'

	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	timeNR, hdot_plus, hdot_cross = np.loadtxt(file_location+datafile_hNRdot, unpack=True)
	timeNR, h_plus, h_cross = np.loadtxt(file_location+datafile_hNR, unpack=True)


	#Making plots
	legend_size = 1
	fig = plt.figure()
	plt.suptitle('Memory Calculation for aligned Sz spin ='+str(mass_ratio_vec[i]) ,fontsize = 15)
	fontP = FontProperties()
	fontP.set_size('10.')


	legend = plt.legend(loc='best',prop={'size':legend_size})
	plt.subplot(2,2,1)
	plt.plot(timeNR, hdot_plus, 'r', label=r'$hdot_{+}$')
	plt.plot(timeNR, hdot_cross, 'k--', label=r'$hdot_{\times}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$hdot_{+,\times}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	plt.subplot(2,1,2)
	plt.plot(timeNR, h_plus, 'r', label=r'$h_{+}$')
	plt.plot(timeNR, h_cross, 'k--', label=r'$h_{\times}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$h_{+,\times}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	#legend = plt.legend(loc='best',prop={'size':legend_size})
	plt.subplot(2,2,2)
	plt.plot(timeNR, hmem, 'r', label=r'$h_{mem}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$h_{mem}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_alignedSpin/'+filename+'.png')
	i+=1
	plt.close()



i=0
#Making plots
legend_size = 2
fig = plt.figure()
plt.title('Memory Calculation for aligned Sz spin ' ,fontsize = 15)
fontP = FontProperties()
fontP.set_size('10.')

legend = plt.legend(loc='best',prop={'size':legend_size})

for filename in filename_vec:
	
	datafile_hMemNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hMemNR.dat'
	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	plt.plot(timeNR, hmem, label=filename_vec[i])
	i+=1

plt.grid()
#plt.ylim(0,0.0007)
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{mem}$')
plt.legend(loc=2)
fontP.set_size('10.')
plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_alignedSpin/'+filename+'.pdf')
plt.show()
	
