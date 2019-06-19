import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
#import data
filename_vec=['0p20', '0p60', '0p80','0p99']
Spin_vec = [0.20, 0.60, 0.80, 0.99]


#NonSpinning data

file_nospin ='/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/Memory_data/rMPsi4_noSpin_q1dataClean_hMemNR.dat'
#import data
	
timeNR_nospin, hmem_nospin, h_mem_plus_nospin = np.loadtxt(file_nospin, unpack=True)
#Normalize to stich
hmem_nospin*=17.0


filename_vec=['m0p20', 'm0p43', 'm0p60','m0p90']
Spin_vec = [-0.20, -0.43, -0.60, -0.90]


                                         
# defining intial time for given solar mass for a given fixed time interval of two weeks
def time_initial(x):
	return 25.0 - 2.6*pow(10, 11-x)


# Function that compute memory growth during two weeks interval of time for a given mass

def Memory_growth_in_two_weeks(log_SolarMass, Spin):
	t_initial =  time_initial(log_SolarMass)
	
	if Spin < 0:
		file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
		filename = 'm0p'+str(Spin)[3:5]
		
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	
	if Spin == 0:
		file_location ='/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/Memory_data/rMPsi4_noSpin_q1dataClean_hMemNR.dat'

		filename = 'm0p'+str(Spin)[2:4]
		
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	
	if Spin > 0:
		file_location ='/home/ashok/gravitational_wave_memory_project/data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
		filename = '0p'+str(Spin)[2:4]
	
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)

	idx = find_nearest_idx(timeNR, t_initial)
	
	hmem_diff = abs(hmem[find_nearest_idx(timeNR, 25)]-hmem[idx])

	return hmem_diff

		
for i in range(8,11):
	print i, Memory_growth_in_two_weeks(i, -0.204), Memory_growth_in_two_weeks(i, 0.604), Memory_growth_in_two_weeks(i, 0.804)
        
