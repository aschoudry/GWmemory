import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties
from matplotlib import colors, cm
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm
import plotsettings
from matplotlib.font_manager import FontProperties


def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
                                         
# defining intial time for given solar mass for a given fixed time interval of two weeks
def time_initial(x, time_final):
	return time_final - 2.6*pow(10, 11-x)

#function that integrate the memory signal
def s(t, hpmem, mu, phi):
	dt =t[1]-t[0]
	int_hp = np.cumsum(hpmem)*dt
	
	st = (1.0/2.0)*(1.0+mu)*(int_hp*np.cos(2*phi))

	return st


# Function that compute memory growth during two weeks interval of time for a given mass

def Memory_growth_in_two_weeks(log_SolarMass, Spin):
		
	if Spin < 0:
		file_location_hmem ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
		filename = 'm0p'+str(Spin)[3:5]
		
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
	
		datafile_hNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hNR.dat'
		timeNR, hp, hc = np.loadtxt(file_location_hmem+datafile_hNR, unpack=True)
	

		
	if Spin == 0.0:
		file_location_hmem ='/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/Memory_data/'

		filename_hmem = 'rMPsi4_noSpin_q1dataClean_hMemNR.dat'
		filename_h22 = 'rMPsi4_noSpin_q1dataClean_hNR.dat'

		datafile_hMemNR = filename_hmem
		datafile_hNR = filename_h22
		
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
		timeNR, hp, hc = np.loadtxt(file_location_hmem+datafile_hNR, unpack=True)	

	if Spin > 0:
		file_location_hmem ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
		filename = '0p'+str(Spin)[2:4]
	
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)

		datafile_hNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hNR.dat'
		timeNR, hp, hc = np.loadtxt(file_location_hmem+datafile_hNR, unpack=True)

	hmem_diff = np.zeros(len(hmem))
	for i in range(len(hmem)):
		hmem_diff[i] =  (hmem[-i]-hmem[-(i+1)])/abs(hmem[-i])

	
	#time_final=timeNR[np.argmax(hp)+135]
	#time_final=timeNR[np.argmin(np.diff(np.diff(hmem)))]
	time_final = timeNR[-find_nearest_idx(hmem_diff, 10.0**(-5))]

	t_initial =  time_initial(log_SolarMass, time_final)

	idx = find_nearest_idx(timeNR, t_initial)

	hmem_two_weeks = hmem[idx:find_nearest_idx(timeNR, time_final)]
	time_two_weeks = timeNR[idx:find_nearest_idx(timeNR, time_final)]
	

	return timeNR, hmem, time_two_weeks, hmem_two_weeks


# Make heatmap for the two weeks memory growth in SMBHB		
#Spin_array = np.array([-0.901, -0.601, -0.431, -0.201, 0.0 ,0.201, 0.601, 0.801, 0.99])

legend_size = 2
fig = plt.figure()
fontP = FontProperties()
fontP.set_size('20.')

Spin_array = np.array([-0.901,0.0 ,0.99])

Mass_array = np.array([8, 9, 10])


shift = 2500
k=0
for i in Spin_array:
	for j in Mass_array:
			timeNR = Memory_growth_in_two_weeks(j, i)[0]
			hmem = Memory_growth_in_two_weeks(j, i)[1]
			time_two_weeks = Memory_growth_in_two_weeks(j, i)[2]
			hmem_two_weeks = Memory_growth_in_two_weeks(j, i)[3]

			plt.plot(timeNR+k*shift, hmem,'k--')
			plt.plot(time_two_weeks+k*shift, hmem_two_weeks, label='M='+str(round(j,1))+' \t\t S ='+ str(round(i, 2)))
			k+=1		

#plt.xlim(-3500, 4500)
plt.xlabel(r'$t/M$')
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$')
plt.legend(loc=2)
fontP.set_size('12.')
plt.show()

plt.plot(np.diff(np.diff(hmem)))
plt.show()


