import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties


#Load saved memory data from mathematica in geometric units

file_location_hmemGeoUnits = "/home/aschoudhary/gravitational_wave_memory_project/data/dataForPaper/MemoryGrowthInDaysIngral/"
		

# Load two weeks data for different mass values

filename_Spin0p99_10Pow8MSun = 'MemoryGrowthIngrlSpin0p99InDays14LogMass8p.txt'
filename_Spin0p0_10Pow8MSun = 'MemoryGrowthIngrlSpin0p0InDays14LogMass8p.txt'
filename_Spinm0p94_10Pow8MSun = 'MemoryGrowthIngrlSpinm0p94InDays14LogMass8p.txt'

filename_Spin0p99_10Pow10MSun = 'MemoryGrowthIngrlSpin0p99InDays14LogMass10p.txt'
filename_Spin0p0_10Pow10MSun = 'MemoryGrowthIngrlSpin0p0InDays14LogMass10p.txt'
filename_Spinm0p94_10Pow10MSun = 'MemoryGrowthIngrlSpinm0p94InDays14LogMass10p.txt'

time_Spin0p99_10Pow8MSun, hmem_Spin0p99_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow8MSun, unpack=True)
time_Spin0p0_10Pow8MSun, hmem_Spin0p0_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow8MSun, unpack=True)
time_Spinm0p94_10Pow8MSun, hmem_Spinm0p94_10Pow8MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow8MSun, unpack=True)

time_Spin0p99_10Pow10MSun, hmem_Spin0p99_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p99_10Pow10MSun, unpack=True)
time_Spin0p0_10Pow10MSun, hmem_Spin0p0_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spin0p0_10Pow10MSun, unpack=True)
time_Spinm0p94_10Pow10MSun, hmem_Spinm0p94_10Pow10MSun = np.loadtxt(file_location_hmemGeoUnits+filename_Spinm0p94_10Pow10MSun, unpack=True)

### Making Plots ########
legend_size = 5
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')
plt.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.95, wspace= 0.2, hspace = 0.2)

plt.subplot(2,2,1)

plt.plot(time_Spin0p99_10Pow8MSun, hmem_Spin0p99_10Pow8MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.plot(time_Spin0p0_10Pow8MSun, hmem_Spin0p0_10Pow8MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$')
plt.plot(time_Spinm0p94_10Pow8MSun, hmem_Spinm0p94_10Pow8MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$')

plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
plt.xticks([0, 5, 10, 15])
#plt.yticks([2.5, 5])
#plt.ylim(0, 0.1)
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()

#plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDays10pow8MSun.pdf')
#plt.show()

##Calculate the quadratic fit and Subtract

res0p99_10pow8= hmem_Spin0p99_10Pow8MSun 
res0p0_10pow8= hmem_Spin0p0_10Pow8MSun 
resm0p94_10pow8= hmem_Spinm0p94_10Pow8MSun 

res0p99_10pow8_quadfit = np.polyfit(time_Spin0p99_10Pow8MSun, res0p99_10pow8, 2)
res0p99_10pow8_1d = np.poly1d(res0p99_10pow8_quadfit)
quadratic_fit_to_res0p99_10pow8 = res0p99_10pow8_1d(time_Spin0p99_10Pow8MSun)
res0p99_10pow8_quadSubtract = res0p99_10pow8-quadratic_fit_to_res0p99_10pow8


#fig= plt.figure()
#fontP = FontProperties()

plt.subplot(2,2,2)

plt.plot(time_Spin0p99_10Pow8MSun, hmem_Spin0p99_10Pow8MSun, 'b', label='Residual')
plt.plot(time_Spin0p99_10Pow8MSun, quadratic_fit_to_res0p99_10pow8, 'r--', label='Quadratic fit')
plt.plot(time_Spin0p99_10Pow8MSun, res0p99_10pow8_quadSubtract, 'g', label='Residual post quadratric subtraction')

plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
plt.xticks([0, 5, 10, 15])
#plt.ylim(0, 0.1)
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()

#plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksQuadFitInDays10pow8MSun.pdf')

#plt.show()


#fig= plt.figure()
#fontP = FontProperties()

plt.subplot(2,2,3)

plt.plot(time_Spin0p99_10Pow10MSun, hmem_Spin0p99_10Pow10MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.plot(time_Spin0p0_10Pow10MSun, hmem_Spin0p0_10Pow10MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.0$')
plt.plot(time_Spinm0p94_10Pow10MSun, hmem_Spinm0p94_10Pow10MSun, alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =-0.94$')


plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
plt.xticks([0, 5, 10, 15])
#plt.ylim(0, 0.1)
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()


#plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDays10pow10MSun.pdf')
#plt.show()



##Calculate the quadratic fit and Subtract

res0p99_10pow10= hmem_Spin0p99_10Pow10MSun 
res0p0_10pow10= hmem_Spin0p0_10Pow10MSun 
resm0p94_10pow10= hmem_Spinm0p94_10Pow10MSun 

res0p99_10pow10_quadfit = np.polyfit(time_Spin0p99_10Pow10MSun, res0p99_10pow10, 2)
res0p99_10pow10_1d = np.poly1d(res0p99_10pow10_quadfit)
quadratic_fit_to_res0p99_10pow10 = res0p99_10pow10_1d(time_Spin0p99_10Pow10MSun)
res0p99_10pow10_quadSubtract = res0p99_10pow10-quadratic_fit_to_res0p99_10pow10

#fig= plt.figure()
#fontP = FontProperties()

plt.subplot(2,2,4)

plt.plot(time_Spin0p99_10Pow10MSun, hmem_Spin0p99_10Pow10MSun, 'b', label='Residual')
plt.plot(time_Spin0p99_10Pow10MSun, quadratic_fit_to_res0p99_10pow10, 'r--', label='Quadratic fit')
plt.plot(time_Spin0p99_10Pow10MSun, res0p99_10pow10_quadSubtract, 'g', label='Residual post quadratric subtraction')
plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')

plt.xticks([0, 5, 10, 15])
#plt.ylim(0, 0.1)
fontP.set_size('5.0')
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksQuadFitInDays.pdf')
plt.show()

