import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties
from matplotlib import colors, cm
from matplotlib.ticker import LogLocator
from matplotlib.colors import LogNorm

numer_of_observation_days=14.0	

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
        #Increase the resolution which works for 10^12 Msun data
        #'''

        if log_SolarMass > 10.0:

            dt_NR_intrp=(timeNR[2]-timeNR[1])/10.0
            timeNR_intrp = np.arange(timeNR[0], timeNR[-1], dt_NR_intrp)

            hmem_intrp = np.interp(timeNR_intrp, timeNR, hmem)

            timeNR = timeNR_intrp
            hmem = hmem_intrp
        #'''

        #Extend the time period when the memory settlles to make it look like a step function

        timeNR_up=timeNR- timeNR[0]+timeNR[-1]
        timeNR=np.append(timeNR,timeNR_up)

        hmem_up=np.full(len(timeNR_up), max(hmem))
        hmem=np.append(hmem, hmem_up)

        #Length of NR data
        len_NR_data = len(timeNR)

        # Create a comple memory waveform
        idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.8)
        
        #Shifting the Post-Newtonian part to where it blows up
        
        hp_mem_PN = hmem_PN 
        hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
        hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
        time_PN = timePN[:hp_mem_PN_max_idx]
                
        time_PN=time_PN-time_PN[-1]
        time_PN=time_PN+timeNR[idx_NRpeak]

        #Find index of ti in NR on PN 
        idx_ti = find_nearest_idx(time_PN, timeNR[0])
        

        ## Adjusting NR memory to shift
        #hp_mem_PN=hp_mem_PN-hp_mem_PN[0]

        hmem = hmem+hp_mem_PN[idx_ti]
        timeNR=np.append(time_PN[:idx_ti], timeNR)
        hmem = np.append(hp_mem_PN[:idx_ti], hmem)

	#Look where the memory growth looks greatest

        if log_SolarMass > 10.0:
            idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.999)
            time_final_i = timeNR[idx_NRpeak]
        else:	
            time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

        
	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]
        
        '''
        #check if two weeks interval looks correct 
        plt.plot(time_two_weeks_i, hmem_two_weeks_i)
        plt.plot(timeNR, hmem, '--')
        plt.show()
        '''
        if log_SolarMass > 10.0:
            idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.999)
            time_final_ip1 = timeNR[idx_NRpeak-1]
        else:	
            time_final_ip1 = timeNR[-2]

	t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

	idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

	hmem_two_weeks_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
	time_two_weeks_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


	hmem_growth_two_weeks_i = hmem_two_weeks_i[-1]-hmem_two_weeks_i[0]	
	hmem_growth_two_weeks_ip1 = hmem_two_weeks_ip1[-1]-hmem_two_weeks_ip1[0]

	k=0
	
	while hmem_growth_two_weeks_ip1 >= hmem_growth_two_weeks_i:
		k+=1
		time_final_i = timeNR[-(1+k)]

		t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

		idx_i = find_nearest_idx(timeNR, t_initial_i)

		hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
		time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

		time_final_ip1 = timeNR[-(2+k)]

		t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

		idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

		hmem_two_weeks_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
		time_two_weeks_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


		hmem_growth_two_weeks_i = hmem_two_weeks_i[-1]-hmem_two_weeks_i[0]	
		hmem_growth_two_weeks_ip1 = hmem_two_weeks_ip1[-1]-hmem_two_weeks_ip1[0]

		
		
	hmem_two_weeks = hmem_two_weeks_ip1
	time_two_weeks = time_two_weeks_ip1

	hmem_diff = hmem_two_weeks[-1]-hmem_two_weeks[0]

	return hmem_diff


# Make heatmap for the two weeks memory growth in SMBHB		
#Spin_array = np.array([-0.941, -0.801, -0.601, -0.431, -0.201, 0.0 ,0.201, 0.601, 0.801, 0.99])

Spin_array = np.array([0.99, 0.801, 0.601, 0.201, 0.0, -0.201, -0.431, -0.601, -0.801, -0.941])
Mass_array = np.linspace(8.0, 10.5, 2*len(Spin_array))

Memory_growth_M_vs_Spin = np.zeros([len(Spin_array), len(Mass_array)])

#4.4*(10**(Mass_array[j]-23))
for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Memory_growth_M_vs_Spin[i][j]  = np.log10(4.9*(10**(Mass_array[j]-23))*Memory_growth_in_two_weeks(Mass_array[j], Spin_array[i], numer_of_observation_days))

vmn = Memory_growth_M_vs_Spin.min()
vmx = Memory_growth_M_vs_Spin.max()

fig, ax = plt.subplots()
im = plt.imshow(Memory_growth_M_vs_Spin, interpolation='bilinear', cmap=cm.RdYlGn, extent=[8, 12.0, -0.94, 0.99], vmin=vmn, vmax=vmx)

plt.xlabel('Log{(M)')
plt.ylabel('Spin')
#plt.clim(abs(Memory_growth_M_vs_Spin).min(),abs(Memory_growth_M_vs_Spin).max())
plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/Memory_Spin_vs_Mass_v2.pdf")	
#plt.show()

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

	#Look where the memory growth looks greatest

	hmem=17.0*hmem
        #Increase the resolution which works for 10^12 Msun data
        #'''

        if log_SolarMass > 10.0:

            dt_NR_intrp=(timeNR[2]-timeNR[1])/10.0
            timeNR_intrp = np.arange(timeNR[0], timeNR[-1], dt_NR_intrp)

            hmem_intrp = np.interp(timeNR_intrp, timeNR, hmem)

            timeNR = timeNR_intrp
            hmem = hmem_intrp
        #'''


        #Extend the time period when the memory settlles to make it look like a step function
   
        timeNR_up=timeNR- timeNR[0]+timeNR[-1]
        timeNR=np.append(timeNR,timeNR_up)

        hmem_up=np.full(len(timeNR_up), max(hmem))
        hmem=np.append(hmem, hmem_up)

        #Length of NR data
        len_NR_data = len(timeNR)

        # Create a comple memory waveform
        idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.8)
        
        #Shifting the Post-Newtonian part to where it blows up
        
        hp_mem_PN = hmem_PN 
        hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
        hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
        time_PN = timePN[:hp_mem_PN_max_idx]
                
        time_PN=time_PN-time_PN[-1]
        time_PN=time_PN+timeNR[idx_NRpeak]

        #Find index of ti in NR on PN 
        idx_ti = find_nearest_idx(time_PN, timeNR[0])
        

        ## Adjusting NR memory to shift
        #hp_mem_PN=hp_mem_PN-hp_mem_PN[0]

        hmem = hmem+hp_mem_PN[idx_ti]
        timeNR=np.append(time_PN[:idx_ti], timeNR)
        hmem = np.append(hp_mem_PN[:idx_ti], hmem)

	
        if log_SolarMass > 10.0:
            idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.999)
            time_final_i = timeNR[idx_NRpeak]
        else:	
            time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

        if log_SolarMass > 10.0:
            idx_NRpeak = find_nearest_idx(timeNR, max(hmem)*0.999)
            time_final_ip1 = timeNR[idx_NRpeak-1]
        else:	
            time_final_ip1 = timeNR[-2]

	t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

	idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

	hmem_two_weeks_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
	time_two_weeks_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


	hmem_growth_two_weeks_i = hmem_two_weeks_i[-1]-hmem_two_weeks_i[0]	
	hmem_growth_two_weeks_ip1 = hmem_two_weeks_ip1[-1]-hmem_two_weeks_ip1[0]

	k=0
	
	while hmem_growth_two_weeks_ip1 >= hmem_growth_two_weeks_i:
		k+=1
		time_final_i = timeNR[-(1+k)]

		t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

		idx_i = find_nearest_idx(timeNR, t_initial_i)

		hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
		time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

		time_final_ip1 = timeNR[-(2+k)]

		t_initial_ip1 =  time_initial(log_SolarMass, time_final_ip1, numer_of_observation_days)

		idx_ip1 = find_nearest_idx(timeNR, t_initial_ip1)

		hmem_two_weeks_ip1 = hmem[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]
		time_two_weeks_ip1 = timeNR[idx_ip1:find_nearest_idx(timeNR, time_final_ip1)]


		hmem_growth_two_weeks_i = hmem_two_weeks_i[-1]-hmem_two_weeks_i[0]	
		hmem_growth_two_weeks_ip1 = hmem_two_weeks_ip1[-1]-hmem_two_weeks_ip1[0]

		
	timeNR_two_weeks = time_two_weeks_ip1
	hmem_two_weeks =  4.9*(10**(log_SolarMass-23))*hmem_two_weeks_ip1

        
        hmem_two_weeks_up=np.full(len(timeNR_two_weeks), max(hmem_two_weeks))
        hmem_two_weeks=np.append(hmem_two_weeks, hmem_two_weeks_up)

        timeNR_two_weeks_up=timeNR_two_weeks- timeNR_two_weeks[0]+timeNR_two_weeks[-1]
        timeNR_two_weeks=np.append(timeNR_two_weeks, timeNR_two_weeks_up)

        print len(timeNR_two_weeks), len(hmem_two_weeks)
	#Convert to femto seconds
	hmem_two_weeks=pow(10,15)*hmem_two_weeks
	
	#compute the quadratic fit and subtact it from reseduals
	res = s(timeNR_two_weeks, hmem_two_weeks, 0, 0) 
	res_quadfit = np.polyfit(timeNR_two_weeks, res, 2)
	res1d = np.poly1d(res_quadfit)
	quadratic_fit_to_res = res1d(timeNR_two_weeks)
	res_quadSubtract = res-quadratic_fit_to_res

	
	# mean of reseduals
	res_mean = 4.9*(10**(log_SolarMass-23))*np.sqrt(np.mean(res_quadSubtract**2))

	return res_mean
      
# Make heatmap for the two weeks memory growth in SMBHB	
numer_of_observation_days=7.0	

Reseduals_M_vs_Spin = np.zeros([len(Spin_array), len(Mass_array)])

for i in range(len(Spin_array)):
	for j in range(len(Mass_array)):
		Reseduals_M_vs_Spin[i][j]  = np.log10(compute_rms_reseduals(Mass_array[j], Spin_array[i], numer_of_observation_days))

vmn=Reseduals_M_vs_Spin.min()
vmx=Reseduals_M_vs_Spin.max()

fig, ax = plt.subplots()
im = plt.imshow(Reseduals_M_vs_Spin , interpolation='bilinear', cmap=cm.RdYlGn, extent=[8, 12.0, -0.94, 0.99], vmin=vmn, vmax=vmx)

plt.xlabel('Log{(M)')
plt.ylabel('Spin')
plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/MemoryRes_Spin_vs_Mass_v2.pdf")
#plt.show()


print Memory_growth_M_vs_Spin.min(), Memory_growth_M_vs_Spin.max()
print Reseduals_M_vs_Spin.min(), Reseduals_M_vs_Spin.max()	
