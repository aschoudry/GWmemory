import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties
from matplotlib import colors, cm
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm

#Set memory growth rate
grwth = 5.0*pow(10, -5)

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

	hmem=17.0*hmem	

	hmem_diff = np.zeros(len(hmem))
	for i in range(len(hmem)):
		hmem_diff[i] =  (hmem[-i]-hmem[-(i+1)])/abs(hmem[-i])

	
	time_final = timeNR[-find_nearest_idx(hmem_diff, grwth)]	

	t_initial =  time_initial(log_SolarMass, time_final)

	idx = find_nearest_idx(timeNR, t_initial)

#	print timeNR[np.argmax(hp)]	
	hmem_diff = abs(hmem[find_nearest_idx(timeNR, time_final)]-hmem[idx])
	
#	plt.plot(timeNR, hp)
#	plt.show()

	return hmem_diff


# Make heatmap for the two weeks memory growth in SMBHB		
Spin_array = np.array([-0.941, -0.801, -0.601, -0.431, -0.201, 0.0 ,0.201, 0.601, 0.801, 0.99])
Mass_array = np.linspace(8.0, 9.5, len(Spin_array))

Memory_growth_M_vs_Spin = np.zeros([len(Spin_array), len(Spin_array)])

#4.4*(10**(Mass_array[j]-23))
for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Memory_growth_M_vs_Spin[i][j]  = 4.6*(10**(Mass_array[j]-23))*Memory_growth_in_two_weeks(Mass_array[j], Spin_array[i])

fig, ax = plt.subplots()
im = plt.imshow(Memory_growth_M_vs_Spin, interpolation='bilinear', cmap=cm.RdYlGn, origin='lower', extent=[-0.94, 0.99, 8, 9.5], norm=LogNorm(vmin=7.711959855858058e-17, vmax=7.311708502287188e-15))
plt.xlabel('Spin')
plt.ylabel('Log{(M)')
plt.clim(abs(Memory_growth_M_vs_Spin).min(),abs(Memory_growth_M_vs_Spin).max())
plt.colorbar()
fig.tight_layout()
#plt.savefig("/home/ashok/Desktop/gravitational_wave_memory_project/plots/Memory_Spin_vs_Mass.pdf")	
plt.show()


#print abs(Memory_growth_M_vs_Spin).min(),abs(Memory_growth_M_vs_Spin).max()
print Memory_growth_M_vs_Spin.min(), Memory_growth_M_vs_Spin.max()

# Compute RMS reseduals as function of Spin Vs Mass

def compute_rms_reseduals(log_SolarMass, Spin): 

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
	
	hmem=17.0*hmem	

	hmem_diff = np.zeros(len(hmem))
	for i in range(len(hmem)):
		hmem_diff[i] =  (hmem[-i]-hmem[-(i+1)])/abs(hmem[-i])

	
	time_final = timeNR[-find_nearest_idx(hmem_diff, grwth)]

	t_initial =  time_initial(log_SolarMass, time_final)

	idx = find_nearest_idx(timeNR, t_initial)

	hmem_two_weeks =  hmem[idx:find_nearest_idx(timeNR, grwth)]
	timeNR_two_weeks = timeNR[idx:find_nearest_idx(timeNR, grwth)]

	#compute the quadratic fit and subtact it from reseduals
	res = s(timeNR_two_weeks, hmem_two_weeks, 0, 0) 
	res_quadfit = np.polyfit(timeNR_two_weeks, res, 2)
	res1d = np.poly1d(res_quadfit)
	res_quadSubtract = res-res1d(timeNR_two_weeks)

	
	# mean of reseduals
	res_mean = 4.6*(10**(log_SolarMass-23))*np.sqrt(np.mean(res_quadSubtract**2))

#	plt.semilogx(timeNR, h_mem_plus)
#	plt.semilogx(timeNR_two_weeks, hmem_two_weeks, 'k--')
#	plt.plot(timeNR_two_weeks, res_quadSubtract, 'k--')
#	plt.show()
	
	return res_mean
      
# Make heatmap for the two weeks memory growth in SMBHB		

Reseduals_M_vs_Spin = np.zeros([len(Spin_array), len(Spin_array)])

for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Reseduals_M_vs_Spin[i][j]  = compute_rms_reseduals(Mass_array[i], Spin_array[j])

fig, ax = plt.subplots()
im = plt.imshow(Reseduals_M_vs_Spin , interpolation='bilinear', cmap=cm.RdYlGn, origin='lower', extent=[-0.94, 0.99, 8, 9.5], norm=LogNorm(vmin=1.3e-17, vmax=7.0e-16))
plt.xlabel('Spin')
plt.ylabel('Log{(M)')
plt.colorbar()
fig.tight_layout()
#plt.savefig("/home/ashok/Desktop/gravitational_wave_memory_project/plots/MemoryRes_Spin_vs_Mass.pdf")
plt.show()


print Memory_growth_M_vs_Spin.min(), Memory_growth_M_vs_Spin.max()
print Reseduals_M_vs_Spin.min(), Reseduals_M_vs_Spin.max()	
