import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties

def step(x,x0, amp):
	if x<=x0:
		ans=0
	else:
		ans=amp
	return ans 

def s(t, hpmem, mu, phi):
	dt =t[1]-t[0]
	int_hp = np.cumsum(hpmem)*dt
	
	st = (1.0/2.0)*(1.0+mu)*(int_hp*np.cos(2*phi))

	return st


#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SamMassRatio_differentSzSpinOnly/Memory_data/'
#filename_vec=['0p2', '0p43', '0p600', '0p800', '0p95']
filename='0p43'

datafile_hMemNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hMemNR.dat'
timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)

#Reproduce fig1 from Pshirkov paper

time = np.linspace(0,10)
hp = np.array([])
for i in range(len(time)):
	hp=np.append(hp,step(time[i],5,1))

res = s(time, hp, 0, 0) 
res_quadfit = np.polyfit(time, res, 2)
res1d = np.poly1d(res_quadfit)
res_quadSubtract = res-res1d(time)

plt.plot(time, res1d(time), "r--")
plt.plot(time, res, "g--")
plt.plot(time, res_quadSubtract, "b")
plt.show()


####take numerical FFT and plot
#hpf = np.fft.fft(hp)
#freq = np.fft.fftfreq(len(hpf),d=time[2]-time[1])
#plt.plot(freq, hpf)
#plt.show()
#################################################

##construction a step function which looks like hmem from SXS data
maxhmemSXS = max(hmem)
timeleft = np.linspace(timeNR[0], timeNR[-1], len(timeNR))
hmemStepleft = np.array([])
for i in range(len(timeNR)):
	hmemStepleft=np.append(hmemStepleft,step(timeNR[i],timeNR[-100], maxhmemSXS))

timeright = np.linspace(timeNR[-1], -timeNR[0], len(timeNR))
hmemStepright = np.full(len(timeright), maxhmemSXS)

timeStep_tot = np.append(timeleft, timeright)
hmemStep_tot = np.append(hmemStepleft,hmemStepright)
timeNR_tot = np.append(timeNR, timeright)
hmemNR_tot = np.append(hmem,hmemStepright)

plt.plot(timeStep_tot, hmemStep_tot, 'r')
plt.plot(timeNR_tot, hmemNR_tot, 'g')
plt.show()

# Compute resedual for comparision

res1 = s(timeStep_tot, hmemStep_tot, 0, 0) 
res_quadfit1 = np.polyfit(timeStep_tot, res1, 2)
res1d1 = np.poly1d(res_quadfit1)
res_quadSubtract1 = res1-res1d1(timeStep_tot)

res2 = s(timeNR_tot, hmemNR_tot, 0, 0) 
res_quadfit2 = np.polyfit(timeNR_tot, res2, 2)
res1d2 = np.poly1d(res_quadfit2)
res_quadSubtract2 = res2-res1d2(timeNR_tot)


plt.plot(timeStep_tot, res1d1(timeStep_tot), "r--")
plt.plot(timeStep_tot, res1, "g--")
plt.plot(timeStep_tot, res_quadSubtract1, "b")
plt.plot(timeNR_tot, res_quadSubtract2, "r")
plt.show()


#Compute FFT of reseduals

res_quadSubtract1FFT = np.fft.fft(res_quadSubtract1)
freq1 = np.fft.fftfreq(len(res_quadSubtract1),d=time[2]-time[1])

res_quadSubtract2FFT = np.fft.fft(res_quadSubtract2)
freq2 = np.fft.fftfreq(len(res_quadSubtract2),d=time[2]-time[1])

plt.loglog(freq1, abs(res_quadSubtract1FFT)*(time[2]-time[1]), "r")
plt.loglog(freq2, abs(res_quadSubtract2FFT)*(time[2]-time[1]), "g")
plt.show()











