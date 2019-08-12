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
import PostNewtonianMemoryFnc as PNM
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)


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

	#Look where the memory growth looks greatest
	numer_of_observation_days=14.0	
	
	time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

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
														


	return timeNR, hmem, time_two_weeks, hmem_two_weeks


# Make heatmap for the two weeks memory growth in SMBHB		
#Spin_array = np.array([-0.901, -0.601, -0.431, -0.201, 0.0 ,0.201, 0.601, 0.801, 0.99])


#Making plots
legend_size = 2
fig, ax1 = plt.subplots()

fontP = FontProperties()
fontP.set_size('10.')

legend = ax1.legend(loc='best',prop={'size':legend_size})

# PN initial time parameter
t_fPN = 550.0

Spin_array = np.array([-0.941,0.0, 0.99])

Mass_array = np.array([8.0, 9.0,10.0])

shift = 2500
shiftSpin=10000
k=0
l=0
for i in Spin_array:
	for j in Mass_array:
		timeNR = Memory_growth_in_two_weeks(j, i)[0]
		hmem = Memory_growth_in_two_weeks(j, i)[1]
		time_two_weeks = Memory_growth_in_two_weeks(j, i)[2]
		hmem_two_weeks = Memory_growth_in_two_weeks(j, i)[3]
	
		timeNR_Original = timeNR
		timeNR = timeNR+k*shift+l*shiftSpin
		time_two_weeks = time_two_weeks+k*shift+l*shiftSpin

		#Attaching the PostNewtonian part spinning binaries
		t_iPN=timeNR_Original[0]
		dt=timeNR[1]-timeNR[0]
		time_PN=np.arange(t_iPN, t_fPN, dt)
		nsteps=len(time_PN)
		eta=0.25
		R=1.0
		M=1.0
		
		x0=pow(-5.0*M/(256.0*t_iPN*eta), 1.0/4.0)
		chiZa=0.0
		chiZs=i

		#Shifting the Post-Newtonian part to where it blows up
		
		hp_mem_PN = PNM.h_plus_mem(np.pi/2.0, eta, M, R, x0, dt, nsteps, 3, chiZa, chiZs)
		hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
		hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
		time_PN = time_PN[:hp_mem_PN_max_idx]
			
		time_PN=time_PN-time_PN[-1]
		time_PN=time_PN+time_two_weeks[-1]

		

		## Adjusting NR memory to shift
		#hp_mem_PN=hp_mem_PN-hp_mem_PN[0]

		hmem = hmem+hp_mem_PN[0]
		hmem_two_weeks=hmem_two_weeks+hp_mem_PN[0]

		ax1.plot(time_PN, hp_mem_PN, 'r:')		
		ax1.plot(timeNR, hmem,'k--')

		ax1.plot(time_two_weeks, hmem_two_weeks, label='M=10**('+ str(round(j,1))+')\, S ='+str(round(i, 2)))
		k+=1		
	l+=1

plt.plot(time_PN, hp_mem_PN, 'r:', label = "PN memory")		
plt.plot(timeNR, hmem,'k--', label = 'NR memory')

#plt.xlim(-5000, 25000)
plt.ylim(0.0, 0.2)
plt.xlabel(r'$t/M$')
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$')
plt.legend(loc=2)
fontP.set_size('8.')

'''
#Zoon in plot
ax2 = plt.axes([0,0,1,1])
ip = InsetPosition(ax1, [0.40,0.40,0.35,0.35])
ax2.set_axes_locator(ip)
# Mark the region corresponding to the inset axes on ax1 and draw lines
# in grey linking the two axes.
mark_inset(ax1, ax2, loc1=2, loc2=4, fc="none", ec='0.2')
ax2.plot(timeNR, hmem,'k--')
ax2.plot(time_two_weeks, hmem_two_weeks)


plt.xlim(22500, 22530)
plt.ylim(0.0379, 0.068)
'''
plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/MemoryGrowthtwoweeks.pdf")
plt.show()

print max(hmem)


# PN initial time parameter
t_iPN = -1000.0
t_fPN = 830.0

'''
for i in Spin_array:
	j=8
	timeNR = Memory_growth_in_two_weeks(j, i)[0]
	hmem = 17*Memory_growth_in_two_weeks(j, i)[1]
	time_two_weeks = Memory_growth_in_two_weeks(j, i)[2]
	hmem_two_weeks = 17*Memory_growth_in_two_weeks(j, i)[3]

	#Attaching the PostNewtonian part spinning binaries 
	t_iPN = timeNR[0]
	
	dt=timeNR[1]-timeNR[0]
	time_PN=np.arange(t_iPN, t_fPN, dt)
	nsteps=len(time_PN)
	eta=0.25
	R=1.0
	M=1.0
	
	
	x0=pow(-5.0*M/(256.0*t_iPN*eta), 1.0/4.0)
	chiZa=0.0
	chiZs=i

	#Shifting the Post-Newtonian part to where it blows up
	
	hp_mem_PN = PNM.h_plus_mem(np.pi/2.0, eta, M, R, x0, dt, nsteps, 3, chiZa, chiZs)
	hp_mem_PN_max_idx = np.argmax(hp_mem_PN)
	hp_mem_PN = hp_mem_PN[:hp_mem_PN_max_idx]
	time_PN = time_PN[:hp_mem_PN_max_idx]
		
	time_PN=time_PN-time_PN[-1]

	plt.plot(time_PN, hp_mem_PN, label = str(round(i, 2)))			

plt.xlim(-5000.0, 1000.0)
plt.ylim(0.0042, 0.076)
plt.xlabel(r'$t/M$')
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$')
plt.legend(loc=2)
fontP.set_size('12.')
plt.show()

'''



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

	#Look where the memory growth looks greatest

	hmem=17.0*hmem

	numer_of_observation_days=14.0	
	
	time_final_i = timeNR[-1]

	t_initial_i =  time_initial(log_SolarMass, time_final_i, numer_of_observation_days)

	idx_i = find_nearest_idx(timeNR, t_initial_i)

	hmem_two_weeks_i = hmem[idx_i:find_nearest_idx(timeNR, time_final_i)]
	time_two_weeks_i = timeNR[idx_i:find_nearest_idx(timeNR, time_final_i)]

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

		
		
	hmem_two_weeks =  4.6*(10**(log_SolarMass-23))*hmem_two_weeks_ip1
	#Convert to femto seconds
	hmem_two_weeks=pow(10,15)*hmem_two_weeks
	
	timeNR_two_weeks = time_two_weeks_ip1

	#compute the quadratic fit and subtact it from reseduals
	res = s(timeNR_two_weeks, hmem_two_weeks, 0, 0) 
	res_quadfit = np.polyfit(timeNR_two_weeks, res, 2)
	res1d = np.poly1d(res_quadfit)
	quadratic_fit_to_res = res1d(timeNR_two_weeks)
	res_quadSubtract = res-quadratic_fit_to_res

	
	#mean of reseduals
	res_mean = np.sqrt(np.mean(res_quadSubtract**2))

	timeNR_two_weeks = 5.3*pow(10, log_SolarMass-11)*timeNR_two_weeks

	return timeNR_two_weeks, res, quadratic_fit_to_res, res_quadSubtract, res_mean
      
legend_size = 2
fig = plt.figure()
fontP = FontProperties()
fontP.set_size('20.')

shift = 15
shiftSpin=30
k=0
l=0
time_last=0
for i in Spin_array:
	for j in Mass_array:
		timeNR_two_weeks = compute_rms_reseduals(j, i)[0]		
		res_two_weeks = compute_rms_reseduals(j, i)[1]
		res_quad_fit_two_weeks = compute_rms_reseduals(j, i)[2]
		res_quadSubtract_two_weeks = compute_rms_reseduals(j, i)[3]

		timeNR_two_weeks = timeNR_two_weeks - timeNR_two_weeks[0] +  k*shift + l*shiftSpin
		plt.plot(timeNR_two_weeks, res_two_weeks,'k--')
		plt.plot(timeNR_two_weeks, res_quad_fit_two_weeks ,'r--')
		plt.plot(timeNR_two_weeks, res_quadSubtract_two_weeks , label='M=10**(' + str(round(j,1))+')\, S ='+ str(round(i, 2)))
		#time_last = timeNR_two_weeks[-1]+25
		k+=1		
	l+=1

plt.plot(timeNR_two_weeks, res_two_weeks,'k--', label= 'NR memory')
plt.plot(timeNR_two_weeks, res_quad_fit_two_weeks ,'r--', label='Quadratic fit')

plt.xlim(-10, 200)
plt.ylim(-15.0, 200.0)
plt.xlabel(r'$time \, (days)$')
plt.ylabel(r'$Residuals \, (fs)$')
plt.legend(loc=2)
fontP.set_size('12.')
plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/ResedualGrowthtwoweeks.pdf")
plt.show()


