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

file_location_hmemGeoUnits = "/home/aschoudhary/gravitational_wave_memory_project/data/dataForPaper/MemoryGrowthInDays/"
		

# Load two weeks data for different mass values

filename_Spin0p99_10Pow8MSun = 'datahmemGrowthInDaysSpin0p99N14Days8MSun.txt'
filename_Spin0p0_10Pow8MSun = 'datahmemGrowthInDaysSpin0p0N14Days8MSun.txt'
filename_Spinm0p94_10Pow8MSun = 'datahmemGrowthInDaysSpinm0p94N14Days8MSun.txt'

filename_Spin0p99_10Pow9MSun = 'datahmemGrowthInDaysSpin0p99N14Days9MSun.txt'
filename_Spin0p0_10Pow9MSun = 'datahmemGrowthInDaysSpin0p0N14Days9MSun.txt'
filename_Spinm0p94_10Pow9MSun = 'datahmemGrowthInDaysSpinm0p94N14Days9MSun.txt'

filename_Spin0p99_10Pow10MSun = 'datahmemGrowthInDaysSpin0p99N14Days10MSun.txt'
filename_Spin0p0_10Pow10MSun = 'datahmemGrowthInDaysSpin0p0N14Days10MSun.txt'
filename_Spinm0p94_10Pow10MSun = 'datahmemGrowthInDaysSpinm0p94N14Days10MSun.txt'

filename_Spin0p99_10Pow11MSun = 'datahmemGrowthInDaysSpin0p99N14Days11MSun.txt'
filename_Spin0p0_10Pow11MSun = 'datahmemGrowthInDaysSpin0p0N14Days11MSun.txt'
filename_Spinm0p94_10Pow11MSun = 'datahmemGrowthInDaysSpinm0p94N14Days11MSun.txt'

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

#Making plots

fig= plt.figure()
fontP = FontProperties()
fontP.set_size('20.')
legend_size = 2

plt.plot(time_Spin0p99_10Pow8MSun, hmem_Spin0p99_10Pow8MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$', linewidth=2)
plt.plot(time_Spin0p0_10Pow8MSun, hmem_Spin0p0_10Pow8MSun, '--', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$', linewidth=2)
plt.plot(time_Spinm0p94_10Pow8MSun, hmem_Spinm0p94_10Pow8MSun,'-.', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$', linewidth=2)

plt.xlabel(r'$t(days)$', fontsize=18)
plt.ylabel(r'$(1Gpc/D_{L})h^{(mem)}_{+}$', fontsize=18)
#plt.xticks([-3000, -2000, -1000, 0])
#plt.ylim(0, 0.1)
plt.legend()
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/MemoryGrowth2weeksInDays10pow8MSun.pdf')
plt.show()

fig= plt.figure()
fontP = FontProperties()
fontP.set_size('20.')
legend_size = 2

plt.plot(time_Spin0p99_10Pow9MSun, hmem_Spin0p99_10Pow9MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$', linewidth=2)
plt.plot(time_Spin0p0_10Pow9MSun, hmem_Spin0p0_10Pow9MSun, '--', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$', linewidth=2)
plt.plot(time_Spinm0p94_10Pow9MSun, hmem_Spinm0p94_10Pow9MSun, '-.', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$', linewidth=2)


plt.xlabel(r'$t(days)$', fontsize=18)
plt.ylabel(r'$(1Gpc/D_{L})h^{(mem)}_{+}$', fontsize=18)
#plt.xticks([-3000, -2000, -1000, 0])
#plt.ylim(0, 0.1)
plt.legend()
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/MemoryGrowth2weeksInDays10pow9MSun.pdf')
plt.show()


fig= plt.figure()
fontP = FontProperties()
fontP.set_size('20.')
legend_size = 2

plt.plot(time_Spin0p99_10Pow11MSun, hmem_Spin0p99_10Pow11MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$', linewidth=2)
plt.plot(time_Spin0p0_10Pow11MSun, hmem_Spin0p0_10Pow11MSun, '--', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$', linewidth=2)
plt.plot(time_Spinm0p94_10Pow11MSun, hmem_Spinm0p94_10Pow11MSun, '-.', alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$', linewidth=2)


plt.xlabel(r'$t(days)$', fontsize=18)
plt.ylabel(r'$(1Gpc/D_{L})h^{(mem)}_{+}$', fontsize=18)
#plt.xticks([-3000, -2000, -1000, 0])
#plt.ylim(0, 0.1)
plt.legend()
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/MemoryGrowth2weeksInDays10pow11MSun.pdf')
plt.show()



