import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


#Load saved memory data from mathematica in geometric units

file_location_hmemGeoUnits = "/home/aschoudhary/gravitational_wave_memory_project/data/dataForPaper/MemoryGrowthGeoUnit/"
		
filename_Spin0p99 = 'datahmemGrowthGeoUnitSpin0p99.txt'
filename_Spin0p0 = 'datahmemGrowthGeoUnitSpin0p0.txt'
filename_Spinm0p94 = 'datahmemGrowthGeoUnitSpinm0p94.txt'

time_Spin0p99, hmem_Spin0p99 = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99, unpack=True)
time_Spin0p0, hmem_Spin0p0 = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0, unpack=True)
time_Spinm0p94, hmem_Spinm0p94 = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94, unpack=True)

# Load two weeks data for different mass values

filename_Spin0p99_10Pow8MSun = 'datahmemGrowthGeoUnitSpin0p99TwoWeeks10pow8MSun.txt'
filename_Spin0p0_10Pow8MSun = 'datahmemGrowthGeoUnitSpin0p0TwoWeeks10pow8MSun.txt'
filename_Spinm0p94_10Pow8MSun = 'datahmemGrowthGeoUnitSpinm0p94TwoWeeks10pow8MSun.txt'

filename_Spin0p99_10Pow9MSun = 'datahmemGrowthGeoUnitSpin0p99TwoWeeks10pow9MSun.txt'
filename_Spin0p0_10Pow9MSun = 'datahmemGrowthGeoUnitSpin0p0TwoWeeks10pow9MSun.txt'
filename_Spinm0p94_10Pow9MSun = 'datahmemGrowthGeoUnitSpinm0p94TwoWeeks10pow9MSun.txt'

filename_Spin0p99_10Pow10MSun = 'datahmemGrowthGeoUnitSpin0p99TwoWeeks10pow10MSun.txt'
filename_Spin0p0_10Pow10MSun = 'datahmemGrowthGeoUnitSpin0p0TwoWeeks10pow10MSun.txt'
filename_Spinm0p94_10Pow10MSun = 'datahmemGrowthGeoUnitSpinm0p94TwoWeeks10pow10MSun.txt'

filename_Spin0p99_10Pow11MSun = 'datahmemGrowthGeoUnitSpin0p99TwoWeeks10pow11MSun.txt'
filename_Spin0p0_10Pow11MSun = 'datahmemGrowthGeoUnitSpin0p0TwoWeeks10pow11MSun.txt'
filename_Spinm0p94_10Pow11MSun = 'datahmemGrowthGeoUnitSpinm0p94TwoWeeks10pow11MSun.txt'

time_Spin0p99_10Pow8MSun, hmem_Spin0p99_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow8MSun, unpack=True)
time_Spin0p0_10Pow8MSun, hmem_Spin0p0_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow8MSun, unpack=True)
time_Spinm0p94_10Pow8MSun, hmem_Spinm0p94_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow8MSun, unpack=True)

time_Spin0p99_10Pow9MSun, hmem_Spin0p99_10Pow9MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow9MSun, unpack=True)
time_Spin0p0_10Pow9MSun, hmem_Spin0p0_10Pow9MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow9MSun, unpack=True)
time_Spinm0p94_10Pow9MSun, hmem_Spinm0p94_10Pow9MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow9MSun, unpack=True)

time_Spin0p99_10Pow10MSun, hmem_Spin0p99_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow10MSun, unpack=True)
time_Spin0p0_10Pow10MSun, hmem_Spin0p0_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow10MSun, unpack=True)
time_Spinm0p94_10Pow10MSun, hmem_Spinm0p94_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow10MSun, unpack=True)

time_Spin0p99_10Pow11MSun, hmem_Spin0p99_10Pow11MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow11MSun, unpack=True)
time_Spin0p0_10Pow11MSun, hmem_Spin0p0_10Pow11MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow11MSun, unpack=True)
time_Spinm0p94_10Pow11MSun, hmem_Spinm0p94_10Pow11MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow11MSun, unpack=True)


fig= plt.figure()
fontP = FontProperties()

plt.semilogy(time_Spin0p99, hmem_Spin0p99, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
#plt.semilogy(time_Spin0p0, hmem_Spin0p0, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$')
plt.semilogy(time_Spinm0p94, hmem_Spinm0p94, '--', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$')


plt.scatter([time_Spin0p99_10Pow8MSun[0], time_Spin0p99_10Pow8MSun[-1]], [hmem_Spin0p99_10Pow8MSun[0], hmem_Spin0p99_10Pow8MSun[-1]], marker='x', color='red', label=r'$M=10 ^{8}\ M_\odot \ $')
plt.scatter([time_Spin0p99_10Pow9MSun[0], time_Spin0p99_10Pow9MSun[-1]], [hmem_Spin0p99_10Pow9MSun[0], hmem_Spin0p99_10Pow9MSun[-1]], marker='*', color='black', label=r'$M=10 ^{9}\ M_\odot \ $')
plt.scatter([time_Spin0p99_10Pow11MSun[0], time_Spin0p99_10Pow11MSun[-1]], [hmem_Spin0p99_10Pow11MSun[0], hmem_Spin0p99_10Pow11MSun[-1]], marker='.', color='red', label=r'$M=10 ^{11}\ M_\odot \ $')

#plt.scatter([time_Spin0p0_10Pow8MSun[0], time_Spin0p0_10Pow8MSun[-1]], [hmem_Spin0p0_10Pow8MSun[0], hmem_Spin0p0_10Pow8MSun[-1]], marker='x', color='red')
#plt.scatter([time_Spin0p0_10Pow9MSun[0], time_Spin0p0_10Pow9MSun[-1]], [hmem_Spin0p0_10Pow9MSun[0], hmem_Spin0p0_10Pow9MSun[-1]], marker='*', color='black')
#plt.scatter([time_Spin0p0_10Pow10MSun[0], time_Spin0p0_10Pow10MSun[-1]], [hmem_Spin0p0_10Pow10MSun[0], hmem_Spin0p0_10Pow10MSun[-1]], marker='.', color='red')

plt.scatter([time_Spinm0p94_10Pow8MSun[0], time_Spinm0p94_10Pow8MSun[-1]], [hmem_Spinm0p94_10Pow8MSun[0], hmem_Spinm0p94_10Pow8MSun[-1]], marker='x', color='red')
plt.scatter([time_Spinm0p94_10Pow9MSun[0], time_Spinm0p94_10Pow9MSun[-1]], [hmem_Spinm0p94_10Pow9MSun[0], hmem_Spinm0p94_10Pow9MSun[-1]], marker='*', color='black')
plt.scatter([time_Spinm0p94_10Pow11MSun[0], time_Spinm0p94_10Pow11MSun[-1]], [hmem_Spinm0p94_10Pow11MSun[0], hmem_Spinm0p94_10Pow11MSun[-1]], marker='.', color='red')

plt.xlabel(r'$t/M$', fontsize=18)
plt.ylabel(r'$(R/M)\,h^{(mem)}_{+}$', fontsize=18)
plt.xticks([-3000, -2000, -1000, 0])
plt.ylim(0, 0.1)
plt.legend()
fig.tight_layout()

#plt.semilogy(time_Spin0p99_10Pow11MSun, hmem_Spin0p99_10Pow11MSun)
#plt.semilogy(time_Spin0p0_10Pow11MSun, hmem_Spin0p0_10Pow11MSun)
#plt.semilogy(time_Spinm0p94_10Pow11MSun, hmem_Spinm0p94_10Pow11MSun)
plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/MemoryGrowth2weeksGeoUnits.pdf')
plt.show()



