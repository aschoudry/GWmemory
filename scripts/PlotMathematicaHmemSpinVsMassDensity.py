import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.signal import savgol_filter as sv

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

# define a funcion that computes the residuals

def Compute_residuals(Spin, Log_Mass, days):
    time = Load_data(Spin, Log_Mass, days)[0]
    res = Load_data(Spin, Log_Mass, days)[1]
    res = sv(res, 51, 3)

    res_quadfit = np.polyfit(time, res, 2)

    res_1d = np.poly1d(res_quadfit)

    quadratic_fit_to_res = res_1d(time)

    res_quadSubtract = res-quadratic_fit_to_res

    # mean of reseduals
    res_mean = np.sqrt(np.mean(res_quadSubtract**2))


    return time, res, res_quadSubtract, res_mean, quadratic_fit_to_res


def Residual_growth14Days(Spin, Log_SolarMass):
    return Load_data(Spin, Log_SolarMass, 14)[1][-1]

### Making Plots ########
legend_size = 5
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')
#plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.95, wspace= 0.2, hspace = 0.1)

#plt.subplot(2,2,1)

plt.loglog(Load_data(0.99, 8.0, 14)[0], Load_data(0.99, 8.0, 14)[1], 'c--', alpha=0.7)
plt.loglog(Load_data(0.0, 8.0, 14)[0], Load_data(0.0, 8.0, 14)[1], 'c',alpha=0.7, label=r'$M=10 ^{8}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 8.0, 14)[0], Load_data(-0.94, 8.0, 14)[1], 'c-.', alpha=0.7)

plt.loglog(Load_data(0.99, 9.090910, 14)[0], Load_data(0.99, 9.090910, 14)[1], 'm--', alpha=0.7)
plt.loglog(Load_data(0.0, 9.090910, 14)[0], Load_data(0.0, 9.090910, 14)[1],'m' ,alpha=0.7, label=r'$M=10 ^{9}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 9.090910, 14)[0], Load_data(-0.94, 9.090910, 14)[1], 'm-.', alpha=0.7)


plt.loglog(Load_data(0.99, 11.0909, 14)[0], Load_data(0.99, 11.0909, 14)[1], 'k--', alpha=0.7)
plt.loglog(Load_data(0.0, 11.0909, 14)[0], Load_data(0.0, 11.0909, 14)[1],'k' ,alpha=0.7, label=r'$M=10 ^{11}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 11.0909, 14)[0], Load_data(-0.94, 11.0909, 14)[1], 'k-.', alpha=0.7)



plt.xlabel('t(days)', fontsize=18)
plt.ylabel('pre-fit Residuals', fontsize=18)
#plt.xticks([0, 5, 10, 15])
plt.yticks([pow(10, -24), pow(10, -18), pow(10, -12)])
#plt.ylim(pow(10, -24), pow(10, -8))
plt.legend(loc='best', prop={'size':15})
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDaysPreFit.pdf')
plt.show()

fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')

#plt.subplot(2,2,2)
plt.loglog(Load_data(0.99, 8.0, 1825)[0], Load_data(0.99, 8.0, 1825)[1], 'c--', alpha=0.7)
plt.loglog(Load_data(0.0, 8.0, 1825)[0], Load_data(0.0, 8.0, 1825)[1],'c', alpha=0.7, label=r'$M=10 ^{8}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 8.0, 1825)[0], Load_data(-0.94, 8.0, 1825)[1], 'c-.', alpha=0.7)

plt.loglog(Load_data(0.99, 9.090910, 1825)[0], Load_data(0.99, 9.090910, 1825)[1], 'm--', alpha=0.7)
plt.loglog(Load_data(0.0, 9.090910, 1825)[0], Load_data(0.0, 9.090910, 1825)[1],'m' ,alpha=0.7, label=r'$M=10 ^{9}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 9.090910, 1825)[0], Load_data(-0.94, 9.090910, 1825)[1], '-.', alpha=0.7)


plt.loglog(Load_data(0.99, 11.0909, 1825)[0], Load_data(0.99, 11.0909, 1825)[1], 'k--', alpha=0.7)
plt.loglog(Load_data(0.0, 11.0909, 1825)[0], Load_data(0.0, 11.0909, 1825)[1], 'k' ,alpha=0.7, label=r'$M=10 ^{11}\ M_\odot \ $')
plt.loglog(Load_data(-0.94, 11.0909, 1825)[0], Load_data(-0.94, 11.0909, 1825)[1], 'k-.', alpha=0.7)



plt.xlabel('t(days)', fontsize=18)
plt.ylabel('post-fit Residuals', fontsize=18)
#plt.xticks([0, 5, 10, 15])
plt.yticks([pow(10, -19), pow(10, -16), pow(10, -13), pow(10, -10)])
#plt.ylim(pow(10, -24), pow(10, -8))
plt.legend(loc='best', prop={'size':15})
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth5YearsDiffSpinsInDaysPrefit.pdf')
plt.show()


#plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDays10pow8MSun.pdf')
#plt.show()
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')


#plt.subplot(2,2,3)

#plt.loglog(Compute_residuals(0.99, 8.0, 14)[0], abs(Compute_residuals(0.99, 8.0, 14)[2]), 'c--')
plt.loglog(Compute_residuals(0.0, 8.0, 14)[0], abs(Compute_residuals(-0.94, 8.0, 14)[2]),'c', label=r'$M=10 ^{8}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 8.0, 14)[0], abs(Compute_residuals(0.0, 8.0, 14)[2]), 'c-.')

#plt.loglog(Compute_residuals(0.99, 9.090910, 14)[0], abs(Compute_residuals(0.99, 9.090910, 14)[2]), 'm--')
plt.loglog(Compute_residuals(0.0, 9.090910, 14)[0], abs(Compute_residuals(-0.94, 9.090910, 14)[2]),'m', label=r'$M=10 ^{9}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 9.090910, 14)[0], abs(Compute_residuals(0.0, 9.090910, 14)[2]), 'm-.')

#plt.loglog(Compute_residuals(0.99, 11.0909, 14)[0], abs(Compute_residuals(0.99, 11.0909, 14)[2]), 'k--')
plt.loglog(Compute_residuals(0.0, 11.0909, 14)[0], abs(Compute_residuals(-0.94, 11.0909, 14)[2]),'k', label=r'$M=10 ^{11}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 11.0909, 14)[0], abs(Compute_residuals(0.0, 11.0909, 14)[2]), 'k-.')

plt.loglog(Compute_residuals(0.0, 12.0, 14)[0], abs(Compute_residuals(-0.94, 12.0, 14)[2]),'r', label=r'$M=10 ^{12}\ M_\odot \ $')


plt.xlabel('t(days)', fontsize=18)
plt.ylabel('pre-fit Residuals', fontsize=18)
#plt.xticks([0, 5, 10, 15])
plt.yticks([pow(10, -22), pow(10, -19), pow(10, -16)])
plt.legend(loc='best', prop={'size':15})
fig.tight_layout()

#fig= plt.figure()
#fontP = FontProperties()
plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth2weeksDiffSpinsInDaysPostfit.pdf')
plt.show()


#plt.subplot(2,2,4)
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')

#plt.loglog(Compute_residuals(0.99, 8.0, 1825)[0], abs(Compute_residuals(0.99, 8.0, 1825)[2]), 'c--')
plt.loglog(Compute_residuals(0.0, 8.0, 1825)[0], abs(Compute_residuals(-0.94, 8.0, 1825)[2]), 'c', label=r'$M=10 ^{8}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 8.0, 1825)[0], abs(Compute_residuals(0.0, 8.0, 1825)[2]), 'c-.')

#plt.loglog(Compute_residuals(0.99, 9.090910, 1825)[0], abs(Compute_residuals(0.99, 9.090910, 1825)[2]), 'm--')
plt.loglog(Compute_residuals(0.0, 9.090910, 1825)[0], abs(Compute_residuals(-0.94, 9.0909100, 1825)[2]),'m', label=r'$M=10 ^{9}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 9.090910, 1825)[0], abs(Compute_residuals(0.0, 9.090910, 1825)[2]), 'm-.')

#plt.loglog(Compute_residuals(0.99, 11.0909, 1825)[0], abs(Compute_residuals(0.99,11.0909, 1825)[2]), 'k--')
plt.loglog(Compute_residuals(0.0, 11.0909, 1825)[0], abs(Compute_residuals(-0.94, 11.0909, 1825)[2]),'k', label=r'$M=10 ^{11}\ M_\odot \ $')
#plt.loglog(Compute_residuals(-0.94, 11.0909, 1825)[0], abs(Compute_residuals(0.0, 11.0909, 1825)[2]), 'k-.')

plt.loglog(Compute_residuals(0.0, 12.0, 1825)[0], abs(Compute_residuals(-0.94, 12.0, 1825)[2]),'r', label=r'$M=10 ^{12}\ M_\odot \ $')


plt.xlabel('t(days)', fontsize=18)
plt.ylabel('post-fit Residuals', fontsize=18)
#plt.xticks([0, 5, 10, 15])
plt.yticks([pow(10, -19), pow(10, -16), pow(10, -13), pow(10, -10)])
plt.legend(loc='best', prop={'size':15})
fig.tight_layout()

plt.savefig('/home/aschoudhary/gravitational_wave_memory_project/plots/PlotfromMathematicaData/ResidualGrowth5YearsDiffSpinsInDaysPostfit.pdf')
plt.show()

#Make desity plot for memory growth


#SpinVec = np.array([-0.94, -0.90, -0.80, -0.60, -0.43, -0.2, 0.0, 0.2, 0.6, 0.8, 0.99])
SpinVec = np.array([0.99, 0.8, 0.6, 0.2, 0.0, -0.2, -0.43, -0.60, -0.80, -0.90, -0.94])
Log10PowMSunVec = np.array([8.0, 8.18182, 8.36364, 8.54545, 8.72727, 8.90909, 9.09091, 9.27273, 9.45455, 9.63636, 9.81818, 10.0, 10.1818, 10.3636, 10.5455, 10.7273, 10.9091,11.0909])

#, 10.7273, 10.9091, 11.0909, 11.2727, 11.4545, 11.6364, 11.8182, 12.0
'''
# Make heatmap for the two weeks memory growth in SMBHB	
Inthmem_M_vs_Spin14days = np.zeros([len(SpinVec), len(Log10PowMSunVec)])
Inthmem_M_vs_Spin1825days = np.zeros([len(SpinVec), len(Log10PowMSunVec)])

for i in range(len(SpinVec)):
	for j in range(len(Log10PowMSunVec)):
		Reseduals_M_vs_Spin14days[i][j]  = np.log10(Compute_residuals(SpinVec[i], Log10PowMSunVec[j], 14)[3])
                Reseduals_M_vs_Spin1825days[i][j]  = np.log10(Compute_residuals(SpinVec[i], Log10PowMSunVec[j], 1825)[3])


vmn14days=Reseduals_M_vs_Spin14days.min()
vmx14days=Reseduals_M_vs_Spin14days.max()
vmn1825days=Reseduals_M_vs_Spin1825days.min()
vmx1825days=Reseduals_M_vs_Spin1825days.max()
'''


Reseduals_M_vs_Spin14days = np.zeros([len(SpinVec), len(Log10PowMSunVec)])
Reseduals_M_vs_Spin1825days = np.zeros([len(SpinVec), len(Log10PowMSunVec)])

for i in range(len(SpinVec)):
	for j in range(len(Log10PowMSunVec)):
		Reseduals_M_vs_Spin14days[i][j]  = np.log10(Compute_residuals(SpinVec[i], Log10PowMSunVec[j], 14)[3])
                Reseduals_M_vs_Spin1825days[i][j]  = np.log10(Compute_residuals(SpinVec[i], Log10PowMSunVec[j], 1825)[3])


vmn14days=Reseduals_M_vs_Spin14days.min()
vmx14days=Reseduals_M_vs_Spin14days.max()
vmn1825days=Reseduals_M_vs_Spin1825days.min()
vmx1825days=Reseduals_M_vs_Spin1825days.max()

fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')

#fig, ax = plt.subplots()
#plt.subplot(1,2,1)
im = plt.imshow(Reseduals_M_vs_Spin14days , interpolation='quadric', cmap=cm.plasma, extent=[8, 11.0, -0.94, 0.99], vmin=vmn14days, vmax=vmx14days)

plt.xlabel(r'$\log10{(M)}$', fontsize=18)
plt.ylabel(r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N$', fontsize=18)
plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/PlotfromMathematicaData/MemoryResSpinvsMass14Days.pdf")

plt.show()
fig= plt.figure()
fontP = FontProperties()
fontP.set_size('5.0')

#plt.subplot(1,2,2)
im = plt.imshow(Reseduals_M_vs_Spin1825days , interpolation='bicubic', cmap=cm.plasma, extent=[8, 11.0, -0.94, 0.99], vmin=vmn1825days, vmax=vmx1825days)


plt.xlabel(r'$\log10{(M)}$', fontsize=18)
plt.ylabel(r'$\mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N$', fontsize=18)
plt.colorbar()
fig.tight_layout()
plt.savefig("../plots/PlotfromMathematicaData/MemoryResSpinvsMass1825Days.pdf")
plt.show()


#Check why there are contours for high mass, high spin case
'''
file_location_forResiduals14Days = "../plots/PlotfromMathematicaData/ResidualGrowthPlots14Days/"
file_location_forResiduals5Years = "../plots/PlotfromMathematicaData/ResidualGrowthPlots5Years/"

SpinVec_v2 = np.array([0.99, 0.8, 0.6, 0.2, 0.0, -0.2, -0.43, -0.60, -0.80, -0.90, -0.94])

Log10PowMSunVec_v2 = np.array([8.0, 8.18182, 8.36364, 8.54545, 8.72727, 8.90909, 9.09091, 9.27273, 9.45455, 9.63636, 9.81818, 10.0, 10.1818, 10.3636, 10.5455, 10.7273, 10.9091,11.0909])


for i in range(len(SpinVec_v2)):
	for j in range(len(Log10PowMSunVec_v2)):
		time  = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 14)[0]
                res  = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 14)[1]
                res_quadSubtract = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 14)[2]
                quadSubtract = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 14)[4]

                plt.title("$M=10 ^{"+str(Log10PowMSunVec_v2[j])+"}\ M_\odot \ \, \mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N = \, $" +str(SpinVec_v2[i]))
#                plt.plot(time, res, 'c', label="Res pre-fit")
                plt.plot(time, res_quadSubtract, 'k--', label="Res post-fit")
#                plt.plot(time, quadSubtract, 'r-.', label="Quad fit")
                plt.xlabel('$t$(days)')
                plt.ylabel(r'$Residuals$')
                plt.legend()
                fig.tight_layout()
                plt.savefig(file_location_forResiduals14Days+"ResM"+str(Log10PowMSunVec_v2[j]).replace('.', 'p')+"Spin"+str(SpinVec_v2[i]).replace('.', 'p')+".png")
                plt.close()
                #plt.show()


for i in range(len(SpinVec_v2)):
	for j in range(len(Log10PowMSunVec_v2)):
		time  = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 1825)[0]
                res  = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 1825)[1]
                res_quadSubtract = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 1825)[2]
                quadSubtract = Compute_residuals(SpinVec_v2[i], Log10PowMSunVec_v2[j], 1825)[4]

                plt.title("$M=10 ^{"+str(Log10PowMSunVec_v2[j])+"}\ M_\odot \ \, \mathbf{\chi}_{s} \cdot \hat{\mathbf{L}}_N = \, $" +str(SpinVec_v2[i]))
#                plt.plot(time, res, 'c', label="Res pre-fit")
                plt.plot(time, res_quadSubtract, 'k--', label="Res post-fit")
#                plt.plot(time, quadSubtract, 'r-.', label="Quad fit")
                plt.xlabel('$t$(days)')
                plt.ylabel(r'$Residuals$')
                plt.legend()
                fig.tight_layout()
                plt.savefig(file_location_forResiduals5Years+"ResM"+str(Log10PowMSunVec_v2[j]).replace('.', 'p')+"Spin"+str(SpinVec_v2[i]).replace('.', 'p')+".png")
                plt.close()
                #plt.show()
'''
