import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import plotsettings
from matplotlib.font_manager import FontProperties

file_location_hmem_Ingrl ='../data/dataForPaper/MemoryGrowthInDaysIngral/'
datafile_hMem_Ingrl = "MemoryGrowthIngrlSpinm0p9InDays1825LogMass9p81818.txt"

time_InDays, hMem_Ingrl = np.loadtxt(file_location_hmem_Ingrl+datafile_hMem_Ingrl, unpack=True)

plt.plot(time_InDays, hMem_Ingrl)
plt.show()


