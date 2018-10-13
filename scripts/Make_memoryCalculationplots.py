import numpy as np
import matplotlib.pyplot as plt
import plotsettings

#import data
timeNR, hmem, h_mem_plus = np.loadtxt('/home/ashok/gravitational_wave_memory_project/data/hMemNR.dat', unpack=True)
timeNR, hdot_plus, hdot_cross = np.loadtxt('/home/ashok/gravitational_wave_memory_project/data/hdotNR.dat', unpack=True)
timeNR, h_plus, h_cross = np.loadtxt('/home/ashok/gravitational_wave_memory_project/data/hNR.dat', unpack=True)

plt.figure()
plt.plot(timeNR, hdot_plus, 'r', label=r'$hdot_{+}$')
plt.plot(timeNR, hdot_cross, 'k--', label=r'$hdot_{\times}$')
plt.grid()
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{+,\times}$')
plt.legend(loc=3)
plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/hdot_plus_cross.pdf")

plt.figure()
plt.plot(timeNR, h_plus, 'r', label=r'$h_{+}$')
plt.plot(timeNR, h_cross, 'k--', label=r'$h_{\times}$')
plt.grid()
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{+,\times}$')
plt.legend(loc=3)
plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/h_plus_cross.pdf")

plt.figure()
plt.plot(timeNR, hmem, 'r', label=r'$h_{mem}$')
#plt.plot(timeNR, h_mem_plus, 'k', label=r'$h_{mem}+ h_{+}$')
#plt.plot(timeNR, h_plus, 'g', label=r'$h_{+}$')

plt.grid()
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{mem}$')
plt.legend(loc=3)
plt.savefig("/home/ashok/gravitational_wave_memory_project/plots/h_mem.pdf")

