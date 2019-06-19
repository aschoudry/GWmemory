import sys

sys.path.insert(0, '/home/ashok/gravitational_wave_memory_project/src')

import Favata_memoryNonlinearMemoryExpressions as Fvta
import SpinContributionToNonlinearmemory as FvtaS
import matplotlib.pyplot as plt
import numpy as np
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.integrate import odeint
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

###plot lengend setting###
fig = plt.figure()
#plt.title('Memory Calculation for Spinning case' ,fontsize = 15)
fontP = FontProperties()
fontP.set_size('20.')

#legend = plt.legend(loc='best',prop={'size':legend_size})
##############################

t0 = -3000.0
time = np.linspace(t0, 0.0, 10000)
dt=time[1]-time[0]
M=1
R=1
theta = np.pi/2.0
eta=1.0/4.0
PN_order = 2.5
nsteps = len(time)
print len(time)
x0=pow(-5.00*M/(256.0*time[0]*eta), 1.0/4.0)
chiZs=-0.3
chiZa=0.0

x1 = Fvta.Inegrate(x0, dt, nsteps, eta, M, 0)
x2 = Fvta.Inegrate(x0, dt, nsteps, eta, M, 2.5)

x1S = FvtaS.Inegrate(x0, dt, nsteps, eta, M, 2.5, chiZa, chiZs)


# Make plots for NonSpinning case
plt.plot(time, x1, 'r', label="0 PN")
plt.plot(time, x2, 'g', label= "2.5 PN")
plt.xlabel("time")
plt.ylabel("x")
#plt.show()

# Make plots for Spinning case
plt.plot(time, x1S, 'y--', label="2.5 PN Spin")
plt.legend(loc=1)
plt.show()


# Non Spinning case
hp0PN = Fvta.hx_plus_mem(theta, eta, M, R, x1, dt, 0)*np.exp(0j*Fvta.phase(x1, x0, eta, 0))
hp2p5PN = Fvta.hx_plus_mem(theta, eta, M, R, x2, dt, 2.5)*np.exp(0j*Fvta.phase(x2, x0, eta, 2.5))

# Spinning case
hp2p5PNS = FvtaS.hx_plus_mem(theta, eta, M, R, x1S, dt, 2.5, chiZa, chiZs)*np.exp(0j*FvtaS.phase(x1S, x0, eta, 2.5, chiZa, chiZs))

# Plots for Non spining case
plt.plot(time, abs(hp0PN), label= "0 PN")
plt.plot(time, abs(hp2p5PN), label = " 2.5 PN")

#plt.show()

# Plots for  spining case
plt.plot(time, abs(hp2p5PNS), label = "2.5 PN Spin")
plt.xlabel("time")
plt.ylabel("h")
plt.legend(loc=1)
plt.show()

# Plots for Non spining case
plt.plot(x1, abs(hp0PN), label= "0 PN")
plt.plot(x2, abs(hp2p5PN), label = " 2.5 PN")

#plt.show()

# Plots for  spining case
plt.plot(x1S, abs(hp2p5PNS), label = "2.5 PN Spin")
plt.legend(loc=1)
plt.show()


#take FFT non spinning case
hpf1 = np.fft.fft(hp0PN)
hpf2 = np.fft.fft(hp2p5PN)
freq = np.fft.fftfreq(len(hpf1),d=dt)

#take FFT spinning case
hpf1S = np.fft.fft(hp2p5PNS)
freqS = np.fft.fftfreq(len(hpf1S),d=dt)

# Make frequency domain plot Non Spinning case
plt.loglog(freq, abs(hpf1*dt), 'r', label="0 PN")
plt.loglog(freq, abs(hpf2*dt), 'k', label="2.5 PN")
#plt.show()


# Make frequency domain plot Spinning case 
plt.loglog(freqS, abs(hpf1S*dt), 'y--', label = "2.5 PN Spin")
plt.legend(loc=1)
plt.show()




