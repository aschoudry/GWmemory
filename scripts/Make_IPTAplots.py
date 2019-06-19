import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.integrate import odeint
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/NonSpinning_differentMassRatio/Memory_data/'
#import data
mass_ratio_vec = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9.5]
filename_vec=['q1', 'q1p5', 'q2', 'q2p5','q3','q4','q5','q6', 'q7', 'q8','q9p5']

filename = filename_vec[1]

datafile_hNRdot='rMPsi4_noSpin_'+filename+'dataClean_hdotNR.dat'
datafile_hNR='rMPsi4_noSpin_'+filename+'dataClean_hNR.dat'
datafile_hMemNR='rMPsi4_noSpin_'+filename+'dataClean_hMemNR.dat'

timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
timeNR, hdot_plus, hdot_cross = np.loadtxt(file_location+datafile_hNRdot, unpack=True)
timeNR, h_plus, h_cross = np.loadtxt(file_location+datafile_hNR, unpack=True)

#Making plots
legend_size = 5
fig = plt.figure()
fontP = FontProperties()

plt.plot(timeNR, h_plus)
plt.plot(timeNR, h_plus + (0.0002*(timeNR+1000))**2, 'r--')
plt.xlim(-1000, 200)
plt.ylim(-0.5,0.5)
plt.grid()
plt.xlabel(r'$t/M$')
plt.ylabel(r'$(R/M)\,h}_{+}$')
plt.legend(loc=2)
fontP.set_size('13.')

plt.savefig("/home/ashok/Desktop/IPTA slides/Strain.png")
plt.show()

