import numpy as np
import matplotlib.pyplot as plt

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def Resample_ten_years_memory_data(Spin):
		
    if Spin < 0:
            file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAntialigned_27Dec/Memory_data/'
            filename = 'm0p'+str(Spin)[3:5]
                    
            datafile_PN_mem='PN_data_'+filename+'.txt'
            timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)
            
    if Spin == 0.0:
            file_location_hmem ='../data/NonSpinning_differentMassRatio/Memory_data/'

            datafile_PN_mem='PN_data_0.txt'
            
            timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)

    if Spin > 0:
            file_location_hmem ='../data/SXSdata/Spinning_binary_with_SpinAligned_27Dec/Memory_data/'
            filename = '0p'+str(Spin)[2:4]
    
            datafile_PN_mem='PN_data_'+filename+'.txt'
            timePN, hmem_PN = np.loadtxt(file_location_hmem+datafile_PN_mem, unpack=True)


    return timePN, hmem_PN, file_location_hmem, datafile_PN_mem

Spin=-0.941
#Spin_array = np.array([0.99, 0.801, 0.601, 0.201, 0.0, -0.201, -0.431, -0.601, -0.801, -0.941])


time = Resample_ten_years_memory_data(Spin)[0]
hmem = Resample_ten_years_memory_data(Spin)[1]

t0=time[0]
tf_PN=time[-1]

num_of_subintervl = 40.0
len_each_intrvl = (tf_PN-t0)/num_of_subintervl


dt=time[2]-time[1]

time_intr = np.array([])

ti=t0
tf=t0+len_each_intrvl
i=1

while tf< -9000.0:
    dx_sample = len_each_intrvl/(10*i)
    time_updated=np.arange(ti, tf, dx_sample)
    time_intr = np.append(time_intr, time_updated) 
    ti=tf
    tf=ti+len_each_intrvl
    print i, dx_sample, tf


    i+=1
    

hmem_intrp = np.interp(time_intr, time,hmem)


plt.plot(time, hmem)
plt.plot(time_intr, hmem_intrp, '--')
plt.show()

time_stich_idx = find_nearest_idx(time, time_intr[-1])
time_intr = np.append(time_intr, time[time_stich_idx:])
hmem_intrp = np.append(hmem_intrp, hmem[time_stich_idx:])

plt.plot(time, hmem)
plt.plot(time_intr, hmem_intrp, '--')
plt.show()

f=open(Resample_ten_years_memory_data(Spin)[2]+"Resampled"+Resample_ten_years_memory_data(Spin)[3], "w")

for i in range(len(time_intr)):
    str="%.15e \t %.15e \n" % (time_intr[i], hmem_intrp[i])
    f.write(str)

f.close()

'''
#Resample and save the new data
Spin_array = np.array([0.99, 0.801, 0.601, 0.201, 0.0, -0.201, -0.431, -0.601, -0.801, -0.941])

for k in Spin_array:
    print k
    time = Resample_ten_years_memory_data(k)[0]
    hmem = Resample_ten_years_memory_data(k)[1]

    t0=time[0]
    tf_PN=time[-1]

    num_of_subintervl = 50.0
    len_each_intrvl = (tf_PN-t0)/num_of_subintervl

    dt=time[2]-time[1]

    time_intr = np.array([])

    ti=t0
    tf=t0+len_each_intrvl
    i=1

    while tf< -4000.0:
        dx_sample = len_each_intrvl/(10*i)
        print i, dx_sample, tf
        time_updated=np.arange(ti, tf, dx_sample)
        time_intr = np.append(time_intr, time_updated) 
        ti=tf
        tf=ti+len_each_intrvl

        i+=1
        

    hmem_intrp = np.interp(time_intr, time,hmem)


    plt.plot(time, hmem)
    plt.plot(time_intr, hmem_intrp, '--')
    plt.show()

    time_stich_idx = find_nearest_idx(time, time_intr[-1])
    time_intr = np.append(time_intr, time[time_stich_idx:])
    hmem_intrp = np.append(hmem_intrp, hmem[time_stich_idx:])

    f=open(Resample_ten_years_memory_data(k)[2]+"Resampled"+Resample_ten_years_memory_data(k)[3], "w")

    for i in range(len(time_intr)):
        str="%.15e \t %.15e \n" % (time_intr[i], hmem_intrp[i])
        f.write(str)

    f.close()
'''
