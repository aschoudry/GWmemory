import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties
from matplotlib import colors, cm
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm

numer_of_observation_days=50.0	

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx
                                         
# defining intial time for given solar mass for a given fixed time interval of two weeks
def time_initial(log_solarMass, time_final, numer_of_observation_days):
	return time_final - 1.76*pow(10, 10-log_solarMass)*numer_of_observation_days

#function that integrate the memory signal
def s(t, hpmem, mu, phi):
	dt =t[1]-t[0]
	int_hp = np.cumsum(hpmem)*dt
	
	st = (1.0/2.0)*(1.0+mu)*(int_hp*np.cos(2*phi))

	return st


# Function that compute memory growth during two weeks interval of time for a given mass

def Memory_growth_in_two_weeks(log_SolarMass, Spin, numer_of_observation_days):
		
	if Spin < 0:
		file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
		filename = 'm0p'+str(Spin)[3:5]
		
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
	
		datafile_PN_mem='ResampledPN_data_'+filename+'.txt'
		timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)
	

		
	if Spin == 0.0:
		file_location_hmem ='../data/NonSpinning_differentMassRatio/Memory_data/'

		filename_hmem = 'rMPsi4_noSpin_q1dataClean_hMemNR.dat'

		datafile_hMemNR = filename_hmem
		datafile_PN_mem='ResampledPN_data_0.txt'
                
                timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
                timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)

	if Spin > 0:
		file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
		filename = '0p'+str(Spin)[2:4]
	
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)

		datafile_PN_mem='ResampledPN_data_'+filename+'.txt'
		timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)

	hmem=17.0*hmem
        
        # Do interpolation on NR data to increase the resolution. Necessory because the loaded data could not handle 10^12 solar mass two weeks time interval.
        dt_NR_intrp=(timeNR[2]-timeNR[1])/4.0
        timeNR_intrp = np.arange(timeNR[0], timeNR[-1], dt_NR_intrp)

        hmem_intrp = np.interp(timeNR_intrp, timeNR, hmem)
        timeNR = timeNR_intrp
        hmem = hmem_intrp

        #Create a compleate waveform by stiching PN with NR memory
        #Attaching the PostNewtonian part spinning binaries
    	#Shifting the Post-Newtonian part to where it blows up
		
        #hp_mem_PN = hmem_PN[-2*len(timeNR):]

        hp_mem_PN = hmem_PN
	hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
		
                
        hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
	time_PN = timePN[:hp_mem_PN_max_idx]
			
	time_PN=time_PN-time_PN[-1]
	#time_PN=time_PN+time_two_weeks[-1]

	## Adjusting NR memory to shift
	#hp_mem_PN=hp_mem_PN-hp_mem_PN[0]
        ##Shift timeNR
        timeNR_stich=timeNR-timeNR[0]
        timeNR_stich=timeNR_stich+time_PN[-len(hmem)]
            
	hmem_stich=hmem+hp_mem_PN[-len(timeNR_stich)] 

        hmem_PN_NR_stich=np.append(hp_mem_PN[:-(len(hmem_stich)+1)],hmem_stich)
        time_PN_NR_stich=np.append(time_PN[:-(len(hmem)+1)],timeNR_stich)
        time_PN_NR_stich=time_PN_NR_stich-time_PN_NR_stich[-1]

        hmem = hmem_PN_NR_stich
        timeNR = time_PN_NR_stich

	#Look where the memory growth looks greatest
	#numer_of_observation_days=14.0	
	
	time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_numer_of_observation_days_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_numer_of_observation_days_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

	time_final_ip1 = timeNR[-2]

	t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

	idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

	hmem_numer_of_observation_days_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
	time_numer_of_observation_days_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


	hmem_growth_numer_of_observation_days_i = hmem_numer_of_observation_days_i[-1]-hmem_numer_of_observation_days_i[0]	
	hmem_growth_numer_of_observation_days_ip1 = hmem_numer_of_observation_days_ip1[-1]-hmem_numer_of_observation_days_ip1[0]

	k=0
	
	while hmem_growth_numer_of_observation_days_ip1 >= hmem_growth_numer_of_observation_days_i:
		k+=1
		time_final_i = timeNR[-(1+k)]

		t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

		idx_i = find_nearest_idx(timeNR, t_initial_i)

		hmem_numer_of_observation_days_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
		time_numer_of_observation_days_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

		time_final_ip1 = timeNR[-(2+k)]

		t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

		idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

		hmem_numer_of_observation_days_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
		time_numer_of_observation_days_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


		hmem_growth_numer_of_observation_days_i = hmem_numer_of_observation_days_i[-1]-hmem_numer_of_observation_days_i[0]	
		hmem_growth_numer_of_observation_days_ip1 = hmem_numer_of_observation_days_ip1[-1]-hmem_numer_of_observation_days_ip1[0]

		
		
	hmem_numer_of_observation_days = hmem_numer_of_observation_days_ip1
	time_numer_of_observation_days = time_numer_of_observation_days_ip1

	hmem_diff = hmem_numer_of_observation_days[-1]-hmem_numer_of_observation_days[0]

	return hmem_diff


# Make heatmap for the two weeks memory growth in SMBHB		
Spin_array = np.array([0.99, 0.801, 0.601, 0.201, 0.0, -0.201, -0.431, -0.601, -0.801, -0.941])
Mass_array = np.linspace(8.0, 10.0, 2*len(Spin_array))

Memory_growth_M_vs_Spin = np.zeros([len(Spin_array), len(Mass_array)])

#4.4*(10**(Mass_array[j]-23))
for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Memory_growth_M_vs_Spin[i][j]  = np.log10(4.9*(10**(Mass_array[j]-23))*Memory_growth_in_two_weeks(Mass_array[j], Spin_array[i], numer_of_observation_days))

vmn = Memory_growth_M_vs_Spin.min()
vmx = Memory_growth_M_vs_Spin.max()

fig, ax = plt.subplots()
im = plt.imshow(Memory_growth_M_vs_Spin, interpolation='bilinear', cmap=cm.RdYlGn, extent=[8, 10.0, -0.94, 0.99], vmin=vmn, vmax=vmx)

plt.xlabel('Log{(M)')
plt.ylabel('Spin')

#plt.clim(abs(Memory_growth_M_vs_Spin).min(),abs(Memory_growth_M_vs_Spin).max())
plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/Memory_Spin_vs_Mass_tenyears.pdf")	
plt.show()

# Compute RMS reseduals as function of Spin Vs Mass

def compute_rms_reseduals(log_SolarMass, Spin, numer_of_observation_days): 

	if Spin < 0:
		file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
		filename = 'm0p'+str(Spin)[3:5]
		
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
	
		datafile_PN_mem='ResampledPN_data_'+filename+'.txt'
		timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)
	

		
	if Spin == 0.0:
		file_location_hmem ='../data/NonSpinning_differentMassRatio/Memory_data/'

		filename_hmem = 'rMPsi4_noSpin_q1dataClean_hMemNR.dat'

		datafile_hMemNR = filename_hmem
		datafile_PN_mem='ResampledPN_data_0.txt'
                
                timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)
                timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)

	if Spin > 0:
		file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
		filename = '0p'+str(Spin)[2:4]
	
		datafile_hMemNR='rMPsi4_Sz1_'+filename+'_Sz2_'+filename+'_q1p5dataN4Clean_hMemNR.dat'
		timeNR, hmem, h_mem_plus = np.loadtxt(file_location_hmem+datafile_hMemNR, unpack=True)

		datafile_PN_mem='ResampledPN_data_'+filename+'.txt'
		timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)
	
	hmem=17.0*hmem	
        
 
        # Do interpolation on NR data to increase the resolution. Necessory because the loaded data could not handle 10^12 solar mass two weeks time interval.
        dt_NR=timeNR[2]-timeNR[1]
        timeNR_intrp = np.arange(timeNR[0], timeNR[-1], dt_NR/2.0)

        hmem_intrp = np.interp(timeNR_intrp, timeNR, hmem)
        timeNR = timeNR_intrp
        hmem = hmem_intrp


	#Look where the memory growth looks greatest
        #Create a compleate waveform by stiching PN with NR memory
        #Attaching the PostNewtonian part spinning binaries
    	#Shifting the Post-Newtonian part to where it blows up
		
        #hp_mem_PN = hmem_PN[-2*len(timeNR):]

        hp_mem_PN = hmem_PN
	hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
		
                
        hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
	time_PN = timePN[:hp_mem_PN_max_idx]
			
	time_PN=time_PN-time_PN[-1]
	#time_PN=time_PN+time_two_weeks[-1]

	## Adjusting NR memory to shift
	#hp_mem_PN=hp_mem_PN-hp_mem_PN[0]
        ##Shift timeNR
        timeNR_stich=timeNR-timeNR[0]
        timeNR_stich=timeNR_stich+time_PN[-len(hmem)]
            
	hmem_stich=hmem+hp_mem_PN[-len(timeNR_stich)]
        

        hmem_PN_NR_stich=np.append(hp_mem_PN[:-(len(hmem_stich)+1)],hmem_stich)
        time_PN_NR_stich=np.append(time_PN[:-(len(hmem)+1)],timeNR_stich)
        time_PN_NR_stich=time_PN_NR_stich-time_PN_NR_stich[-1]

        hmem = hmem_PN_NR_stich
        timeNR = time_PN_NR_stich

	##
	time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_numer_of_observation_days_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_numer_of_observation_days_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

	time_final_ip1 = timeNR[-2]

	t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

	idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

	hmem_numer_of_observation_days_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
	time_numer_of_observation_days_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


	hmem_growth_numer_of_observation_days_i = hmem_numer_of_observation_days_i[-1]-hmem_numer_of_observation_days_i[0]	
	hmem_growth_numer_of_observation_days_ip1 = hmem_numer_of_observation_days_ip1[-1]-hmem_numer_of_observation_days_ip1[0]

	k=0
	
	while hmem_growth_numer_of_observation_days_ip1 >= hmem_growth_numer_of_observation_days_i:
		k+=1
		time_final_i = timeNR[-(1+k)]

		t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

		idx_i = find_nearest_idx(timeNR, t_initial_i)

		hmem_numer_of_observation_days_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
		time_numer_of_observation_days_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

		time_final_ip1 = timeNR[-(2+k)]

		t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

		idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

		hmem_numer_of_observation_days_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
		time_numer_of_observation_days_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


		hmem_growth_numer_of_observation_days_i = hmem_numer_of_observation_days_i[-1]-hmem_numer_of_observation_days_i[0]	
		hmem_growth_numer_of_observation_days_ip1 = hmem_numer_of_observation_days_ip1[-1]-hmem_numer_of_observation_days_ip1[0]

		
		
	hmem_numer_of_observation_days = hmem_numer_of_observation_days_ip1
	time_numer_of_observation_days = time_numer_of_observation_days_ip1

	timeNR_numer_of_observation_days=time_numer_of_observation_days*4.9*pow(10.0, -6)*pow(10, log_SolarMass)

	#compute the quadratic fit and subtact it from reseduals
	res = s(timeNR_numer_of_observation_days, hmem_numer_of_observation_days, 0, 0) 
	res_quadfit = np.polyfit(timeNR_numer_of_observation_days, res, 2)
	res1d = np.poly1d(res_quadfit)
	res_quadSubtract = res-res1d(timeNR_numer_of_observation_days)

	
	# mean of reseduals
	res_mean = 4.9*(10**(log_SolarMass-23))*np.sqrt(np.mean(res_quadSubtract**2))

	return res_mean
      
# Make heatmap for the two weeks memory growth in SMBHB		

Reseduals_M_vs_Spin = np.zeros([len(Spin_array), len(Mass_array)])

for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Reseduals_M_vs_Spin[i][j]  = np.log10(compute_rms_reseduals(Mass_array[j], Spin_array[i], numer_of_observation_days))

vmn=Reseduals_M_vs_Spin.min()
vmx=Reseduals_M_vs_Spin.max()

fig, ax = plt.subplots()
im = plt.imshow(Reseduals_M_vs_Spin , interpolation='bilinear', cmap=cm.RdYlGn, extent=[8, 10.0, -0.94, 0.99], vmin=vmn, vmax=vmx)

plt.xlabel('Log{(M)')
plt.ylabel('Spin')

plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/MemoryRes_Spin_vs_Mass_tenyears.pdf")
plt.show()


print Memory_growth_M_vs_Spin.min(), Memory_growth_M_vs_Spin.max()
print Reseduals_M_vs_Spin.min(), Reseduals_M_vs_Spin.max()	
