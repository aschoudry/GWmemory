import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties
import memory_function_from_favataPaper as mf

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/Memory_data/'
#import data
mass_ratio_vec = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9.5]
filename_vec=['q1', 'q1p5', 'q2', 'q2p5','q3','q4','q5','q6', 'q7', 'q8','q9p5']
i=0

for filename in filename_vec:
	
	datafile_hNRdot='rMPsi4_noSpin_'+filename+'dataClean_hdotNR.dat'
	datafile_hNR='rMPsi4_noSpin_'+filename+'dataClean_hNR.dat'
	datafile_hMemNR='rMPsi4_noSpin_'+filename+'dataClean_hMemNR.dat'

	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	timeNR, hdot_plus, hdot_cross = np.loadtxt(file_location+datafile_hNRdot, unpack=True)
	timeNR, h_plus, h_cross = np.loadtxt(file_location+datafile_hNR, unpack=True)


	#Making plots
	legend_size = 1
	fig = plt.figure()
	plt.suptitle('Memory Calculation for q ='+str(mass_ratio_vec[i]) ,fontsize = 15)
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

	plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_nonSpining/'+filename+'.png')
	i+=1
	plt.close()




mass_ratio_vec = [2.0, 3.0, 4.0, 5.0, 7.0]
filename_vec=['q2', 'q3', 'q4','q5','q7']




i=0
#Making plots
legend_size = 2
fig = plt.figure()
plt.title('Memory Calculation for NonSpinning ' ,fontsize = 15)
fontP = FontProperties()
fontP.set_size('10.')

legend = plt.legend(loc='best',prop={'size':legend_size})

for filename in filename_vec:
	
	datafile_hMemNR='rMPsi4_noSpin_'+filename+'dataClean_hMemNR.dat'
	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	#Normalize to stich
	hmem*=17.0
	
	#Stiching Postnewtonian memory
	dt=timeNR[1]-timeNR[0]
	time_PN = np.arange(-9000, -2000, dt)
	
	q=mass_ratio_vec[i]
	eta = q/pow(1.0+q,2)
		
	hp_mem_PN = mf.h_plus_mem(np.pi/2, eta, 1.0, 1.0, time_PN)


	dhp_mem_PN_dt=np.diff(hp_mem_PN)/dt

	idx_cut = 80
	dt=timeNR[1]-timeNR[0]
	Slope_hmem_NR=(hmem[idx_cut]-hmem[idx_cut-1])/dt
	idx=find_nearest_idx(dhp_mem_PN_dt, Slope_hmem_NR)
	hmem=hmem[idx_cut:]
	timeNR=timeNR[idx_cut:]
	
	print Slope_hmem_NR, dhp_mem_PN_dt[idx], len(timeNR)

	hp_mem_PN_cut=hp_mem_PN[:idx]
	time_PN_cut=np.linspace(-9000,timeNR[0] , len(hp_mem_PN_cut))


	hmem_tot=hmem+hp_mem_PN_cut[-1]
	time_tot = np.append(time_PN_cut, timeNR)
	hmem_tot = np.append(hp_mem_PN_cut, hmem_tot)

	plt.plot(time_tot, hmem_tot, label=filename_vec[i])
	plt.plot(timeNR, hmem + hp_mem_PN_cut[-1], 'k--' )
	i+=1

plt.grid()
#plt.ylim(0,0.0002)
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{mem}$')
plt.legend(loc=2)
fontP.set_size('10.')
plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_nonSpining/'+filename+'.pdf')
plt.show()
	
