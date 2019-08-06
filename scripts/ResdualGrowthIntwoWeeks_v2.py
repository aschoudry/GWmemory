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

# Set rate of memory growth
grwth = 5.0*pow(10, -5)


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

	
	time_final = timeNR[-find_nearest_idx(hmem_diff, grwth)]

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

# PN initial time parameter
t_fPN = 580.0

Spin_array = np.array([-0.941, 0.0, 0.99])

Mass_array = np.array([8.0, 9.0, 10.0])

shift = 2500
k=0
for i in Spin_array:
	for j in Mass_array:
		timeNR = Memory_growth_in_two_weeks(j, i)[0]
		hmem = 17*Memory_growth_in_two_weeks(j, i)[1]
		time_two_weeks = Memory_growth_in_two_weeks(j, i)[2]
		hmem_two_weeks = 17*Memory_growth_in_two_weeks(j, i)[3]
	
		timeNR_Original = timeNR
		timeNR = timeNR+k*shift
		time_two_weeks = time_two_weeks+k*shift

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
		time_PN=time_PN+timeNR[-1]

		

		## Adjusting NR memory to shift
		#hmem = hmem + hp_mem_PN[-1]
		#hmem_two_weeks = hmem_two_weeks + hp_mem_PN[-1]
		hp_mem_PN=hp_mem_PN-hp_mem_PN[0]
		
		plt.plot(time_PN, hp_mem_PN, 'r:')		

		plt.plot(timeNR, hmem,'k--')
		plt.plot(time_two_weeks, hmem_two_weeks, label='M='+str(round(j,1))+' \t\t S ='+ str(round(i, 2)))
		k+=1		

plt.ylim(0.0, 0.08)
plt.xlabel(r'$t/M$')
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$')
plt.legend(loc=2)
fontP.set_size('12.')
plt.show()

print max(hmem)


# PN initial time parameter
t_iPN = -1000.0
t_fPN = 830.0


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

	hmem_diff = np.zeros(len(hmem))
	for i in range(len(hmem)):
		hmem_diff[i] =  (hmem[-i]-hmem[-(i+1)])/abs(hmem[-i])

	
	time_final = timeNR[-find_nearest_idx(hmem_diff, grwth)]

	t_initial =  time_initial(log_SolarMass, time_final)

	idx = find_nearest_idx(timeNR, t_initial)

	
	#hmem_two_weeks =  4.4*(10**(log_SolarMass-23))*hmem[idx:find_nearest_idx(timeNR, time_final)]
	hmem_two_weeks =  hmem[idx:find_nearest_idx(timeNR, time_final)]
	#hmem = 4.4*(10**(log_SolarMass-23))*hmem
	timeNR_two_weeks = timeNR[idx:find_nearest_idx(timeNR, time_final)]

	#compute the quadratic fit and subtact it from reseduals
	res = s(timeNR_two_weeks, hmem_two_weeks, 0, 0) 
	res_quadfit = np.polyfit(timeNR_two_weeks, res, 2)
	res1d = np.poly1d(res_quadfit)
	quadratic_fit_to_res = res1d(timeNR_two_weeks)
	res_quadSubtract = res-quadratic_fit_to_res

	
	#mean of reseduals
	res_mean = np.sqrt(np.mean(res_quadSubtract**2))

	timeNR_two_weeks = 4.6*pow(10, log_SolarMass-11)*timeNR_two_weeks

	return timeNR_two_weeks, res, quadratic_fit_to_res, res_quadSubtract, res_mean
      
legend_size = 2
fig = plt.figure()
fontP = FontProperties()
fontP.set_size('20.')

shift = 50
k=0
time_last=0
for i in Spin_array:
	for j in Mass_array:
		timeNR_two_weeks = compute_rms_reseduals(j, i)[0]		
		res_two_weeks = compute_rms_reseduals(j, i)[1]
		res_quad_fit_two_weeks = compute_rms_reseduals(j, i)[2]
		res_quadSubtract_two_weeks = compute_rms_reseduals(j, i)[3]

		timeNR_two_weeks = timeNR_two_weeks - timeNR_two_weeks[0] +  k*shift
		plt.plot(timeNR_two_weeks, res_two_weeks,'k--')
		plt.plot(timeNR_two_weeks, res_quad_fit_two_weeks ,'r--')
		plt.plot(timeNR_two_weeks, res_quadSubtract_two_weeks , label='M='+str(round(j,1))+' \t\t S ='+ str(round(i, 2)))
		#time_last = timeNR_two_weeks[-1]+25

		k+=1		

plt.ylim(-0.02, 0.2)
plt.xlabel(r'$days$')
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$')
plt.legend(loc=2)
fontP.set_size('12.')
plt.show()


