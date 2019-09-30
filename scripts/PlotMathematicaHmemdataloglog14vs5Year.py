import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties


#Load saved memory data from mathematica in geometric units

file_location_hmem = "/home/aschoudhary/gravitational_wave_memory_project/data/dataForPaper/MemoryGrowthInDaysIngral/"
		
# Load data for a given input values

def Load_data(Spin, Log_Mass, days):

    if Log_Mass == 8. or Log_Mass ==9. or  Log_Mass ==10. or Log_Mass ==11. or Log_Mass ==12.:
        Log_Mass = str(Log_Mass).replace('.0', 'p')

    else:
        Log_Mass = str(Log_Mass).replace('.', 'p')

    if Spin < 0.0:
        Spin=str(Spin).replace('.', 'p')
        Spin=str(Spin).replace('-', 'm')

    Spin=str(Spin).replace('.', 'p')

    filename = 'MemoryGrowthIngrlSpin'+Spin+'InDays'+str(days)+'LogMass'+Log_Mass+'.txt'

    time, hmem_Intgrl = np.loadtxt(file_location_hmem+filename, unpack=True)
    
    return time, hmem_Intgrl 

# Load two weeks data for different mass values




### Making Plots ########
legend_size = 5
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')
plt.subplots_adjust(left=0.125, right=0.95, bottom=0.1, top=0.95, wspace= 0.2, hspace = 0.2)

plt.subplot(2,2,1)

plt.loglog(Load_data(0.99, 8.0, 14)[0], Load_data(0.99, 8.0, 14)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 8.0, 14)[0], Load_data(0.0, 8.0, 14)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 8.0, 14)[0], Load_data(-0.94, 8.0, 14)[1], '-.', alpha=0.7)

plt.loglog(Load_data(0.99, 9.090910, 14)[0], Load_data(0.99, 9.090910, 14)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 9.090910, 14)[0], Load_data(0.0, 9.090910, 14)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 9.090910, 14)[0], Load_data(-0.94, 9.090910, 14)[1], '-.', alpha=0.7)


plt.loglog(Load_data(0.99, 10.0, 14)[0], Load_data(0.99, 10.0, 14)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 10.0, 14)[0], Load_data(0.0, 10.0, 14)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 10.0, 14)[0], Load_data(-0.94, 10.0, 14)[1], '-.', alpha=0.7)



plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
#plt.xticks([0, 5, 10, 15])
#plt.yticks([2.5, 5])
plt.ylim(pow(10, -24), pow(10, -8))
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()

plt.subplot(2,2,2)
plt.loglog(Load_data(0.99, 8.0, 1825)[0], Load_data(0.99, 8.0, 1825)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 8.0, 1825)[0], Load_data(0.0, 8.0, 1825)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 8.0, 1825)[0], Load_data(-0.94, 8.0, 1825)[1], '-.', alpha=0.7)

plt.loglog(Load_data(0.99, 9.090910, 1825)[0], Load_data(0.99, 9.090910, 1825)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 9.090910, 1825)[0], Load_data(0.0, 9.090910, 1825)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 9.090910, 1825)[0], Load_data(-0.94, 9.090910, 1825)[1], '-.', alpha=0.7)


plt.loglog(Load_data(0.99, 10.0, 1825)[0], Load_data(0.99, 10.0, 1825)[1], '--', alpha=0.7)
plt.loglog(Load_data(0.0, 10.0, 1825)[0], Load_data(0.0, 10.0, 1825)[1], alpha=0.7, label=r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N =0.99$')
plt.loglog(Load_data(-0.94, 10.0, 1825)[0], Load_data(-0.94, 10.0, 1825)[1], '-.', alpha=0.7)



plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
#plt.xticks([0, 5, 10, 15])
#plt.yticks([2.5, 5])
plt.ylim(pow(10, -24), pow(10, -8))
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()


#plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDays10pow8MSun.pdf')
#plt.show()

##Calculate the quadratic fit and Subtract

res0p99_10pow8 = Load_data(0.99, 8.0, 14)[1]
res0p0_10pow8 = Load_data(0.0, 8.0, 14)[1] 
resm0p94_10pow8 = Load_data(-0.94, 8.0, 14)[1] 
time_10Pow8MSun = Load_data(0.99, 8.0, 14)[0]

res0p99_10pow9 = Load_data(0.99, 9.090910, 14)[1]
res0p0_10pow9 = Load_data(0.0, 9.090910, 14)[1] 
resm0p94_10pow9 = Load_data(-0.94, 9.090910, 14)[1] 
time_10Pow9MSun = Load_data(0.99, 9.090910, 14)[0]

res0p99_10pow8 = Load_data(0.99, 8.0, 14)[1]
res0p0_10pow8 = Load_data(0.0, 8.0, 14)[1] 
resm0p94_10pow8 = Load_data(-0.94, 8.0, 14)[1] 
time_10Pow8MSun = Load_data(0.99, 8.0, 14)[0]

res0p99_10pow9 = Load_data(0.99, 9.090910, 14)[1]
res0p0_10pow9 = Load_data(0.0, 9.090910, 14)[1] 
resm0p94_10pow9 = Load_data(-0.94, 9.090910, 14)[1] 
time_10Pow9MSun = Load_data(0.99, 9.090910, 14)[0]


res0p99_10pow8_quadfit = np.polyfit(time_10Pow8MSun, res0p99_10pow8, 2)
res0p0_10pow8_quadfit = np.polyfit(time_10Pow8MSun, res0p0_10pow8, 2)
resm0p94_10pow8_quadfit = np.polyfit(time_10Pow8MSun, resm0p94_10pow8, 2)

res0p99_10pow9_quadfit = np.polyfit(time_10Pow9MSun, res0p99_10pow9, 2)
res0p0_10pow9_quadfit = np.polyfit(time_10Pow9MSun, res0p0_10pow9, 2)
resm0p94_10pow9_quadfit = np.polyfit(time_10Pow9MSun, resm0p94_10pow9, 2)


res0p99_10pow8_1d = np.poly1d(res0p99_10pow8_quadfit)
res0p0_10pow8_1d = np.poly1d(res0p0_10pow8_quadfit)
resm0p94_10pow8_1d = np.poly1d(resm0p94_10pow8_quadfit)

res0p99_10pow9_1d = np.poly1d(res0p99_10pow9_quadfit)
res0p0_10pow9_1d = np.poly1d(res0p0_10pow9_quadfit)
resm0p94_10pow9_1d = np.poly1d(resm0p94_10pow9_quadfit)


quadratic_fit_to_res0p99_10pow8 = res0p99_10pow8_1d(time_10Pow8MSun)
quadratic_fit_to_res0p0_10pow8 = res0p0_10pow8_1d(time_10Pow8MSun)
quadratic_fit_to_resm0p94_10pow8 = resm0p94_10pow8_1d(time_10Pow8MSun)

quadratic_fit_to_res0p99_10pow9 = res0p99_10pow9_1d(time_10Pow9MSun)
quadratic_fit_to_res0p0_10pow9 = res0p0_10pow9_1d(time_10Pow9MSun)
quadratic_fit_to_resm0p94_10pow9 = resm0p94_10pow9_1d(time_10Pow9MSun)

res0p99_10pow8_quadSubtract = res0p99_10pow8-quadratic_fit_to_res0p99_10pow8
res0p0_10pow8_quadSubtract = res0p0_10pow8-quadratic_fit_to_res0p0_10pow8
resm0p94_10pow8_quadSubtract = resm0p94_10pow8-quadratic_fit_to_resm0p94_10pow8

res0p99_10pow9_quadSubtract = res0p99_10pow9-quadratic_fit_to_res0p99_10pow9
res0p0_10pow9_quadSubtract = res0p0_10pow9-quadratic_fit_to_res0p0_10pow9
resm0p94_10pow9_quadSubtract = resm0p94_10pow9-quadratic_fit_to_resm0p94_10pow9

#fig= plt.figure()
#fontP = FontProperties()

plt.subplot(2,2,3)

plt.loglog(time_10Pow8MSun, abs(res0p99_10pow8_quadSubtract), '--')
plt.loglog(time_10Pow8MSun, abs(res0p0_10pow8_quadSubtract),  label='Residual ')
plt.loglog(time_10Pow8MSun, abs(resm0p94_10pow8_quadSubtract), '-.')

plt.loglog(time_10Pow9MSun, abs(res0p99_10pow9_quadSubtract), '--')
plt.loglog(time_10Pow9MSun, abs(res0p0_10pow9_quadSubtract),  label='Residual ')
plt.loglog(time_10Pow9MSun, abs(resm0p94_10pow9_quadSubtract), '-.')

plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
#plt.xticks([0, 5, 10, 15])
#plt.ylim(0, 0.1)
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()
##Calculate the quadratic fit and Subtract

res0p99_10pow8 = Load_data(0.99, 8.0, 1825)[1]
res0p0_10pow8 = Load_data(0.0, 8.0, 1825)[1] 
resm0p94_10pow8 = Load_data(-0.94, 8.0, 1825)[1] 
time_10Pow8MSun = Load_data(0.99, 8.0, 1825)[0]


res0p99_10pow8_quadfit = np.polyfit(time_10Pow8MSun, res0p99_10pow8, 2)
res0p0_10pow8_quadfit = np.polyfit(time_10Pow8MSun, res0p0_10pow8, 2)
resm0p94_10pow8_quadfit = np.polyfit(time_10Pow8MSun, resm0p94_10pow8, 2)


res0p99_10pow8_1d = np.poly1d(res0p99_10pow8_quadfit)
res0p0_10pow8_1d = np.poly1d(res0p0_10pow8_quadfit)
resm0p94_10pow8_1d = np.poly1d(resm0p94_10pow8_quadfit)



quadratic_fit_to_res0p99_10pow8 = res0p99_10pow8_1d(time_10Pow8MSun)
quadratic_fit_to_res0p0_10pow8 = res0p0_10pow8_1d(time_10Pow8MSun)
quadratic_fit_to_resm0p94_10pow8 = resm0p94_10pow8_1d(time_10Pow8MSun)


res0p99_10pow8_quadSubtract = res0p99_10pow8-quadratic_fit_to_res0p99_10pow8
res0p0_10pow8_quadSubtract = res0p0_10pow8-quadratic_fit_to_res0p0_10pow8
resm0p94_10pow8_quadSubtract = resm0p94_10pow8-quadratic_fit_to_resm0p94_10pow8


#fig= plt.figure()
#fontP = FontProperties()

plt.subplot(2,2,4)

plt.loglog(time_10Pow8MSun, abs(res0p99_10pow8_quadSubtract), '--')
plt.loglog(time_10Pow8MSun, abs(res0p0_10pow8_quadSubtract),  label='Residual ')
plt.loglog(time_10Pow8MSun, abs(resm0p94_10pow8_quadSubtract), '-.')



plt.xlabel(r'$t$(days)')
plt.ylabel(r'$Residuals$')
#plt.xticks([0, 5, 10, 15])
#plt.ylim(0, 0.1)
plt.legend(loc='best', prop={'size':legend_size})
fig.tight_layout()
plt.show()

