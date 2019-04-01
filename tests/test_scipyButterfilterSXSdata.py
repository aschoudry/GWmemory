from scipy.signal import butter, lfilter
import numpy as np
from scipy import integrate
from scipy import signal
from scipy import interpolate
import matplotlib.pyplot as plt
import h5py

# Reading the data file for Phi4
file_name_psi="/home/ashok/gravitational_wave_memory_project/data/rMPsi4_Asymptotic_GeometricUnits_CoM.h5"
f_psi = h5py.File(file_name_psi,'r+')

# Reading the strain data
file_name_h="/home/ashok/gravitational_wave_memory_project/data/rhOverM_Asymptotic_GeometricUnits_CoM.h5"
f_h = h5py.File(file_name_h,'r+')

#reading data for time and l2m2 mode
data_psi = f_psi['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
data_h = f_h['Extrapolated_N4.dir']['Y_l2_m2.dat'][:]
time = np.array([])
d2h22_real_by_dt2 = np.array([])
d2h22_imag_by_dt2 = np.array([])
h22_real_SXS=np.array([])
h22_imag_SXS=np.array([])


for i in range(len(data_psi)):
	time=np.append(time, data_psi[:][i][0])
	d2h22_real_by_dt2=np.append(d2h22_real_by_dt2, data_psi[:][i][1])
	d2h22_imag_by_dt2=np.append(d2h22_imag_by_dt2, data_psi[:][i][2])
	h22_real_SXS=np.append(h22_real_SXS,data_h[:][i][1])
	h22_imag_SXS=np.append(h22_imag_SXS,data_h[:][i][2])

plt.plot(time, d2h22_real_by_dt2)
plt.plot(time, d2h22_imag_by_dt2)
plt.show()


dh22_real_by_dt=integrate.cumtrapz(d2h22_real_by_dt2, time,initial=0)
dh22_imag_by_dt=integrate.cumtrapz(d2h22_imag_by_dt2, time, initial=0)

#integrating the data second time with respect to time

h22_real=integrate.cumtrapz(dh22_real_by_dt, time, initial=0)
h22_imag=integrate.cumtrapz(dh22_imag_by_dt, time, initial=0)

time_intrp=np.linspace(time[0], time[-1], len(time))
h22_real_intrp=interpolate.interp1d(time_intrp, h22_real)
h22_imag_intrp=interpolate.interp1d(time_intrp, h22_imag)

time=time_intrp
h22_real=h22_real_intrp(time)
h22_imag=h22_imag_intrp(time)

dt=time[2]-time[1]

def butter_bandpass(lowcut, highcut, fs, order=3):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, lowcut, highcut, fs, order=3):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = lfilter(b, a, data)
    return y


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    from scipy.signal import freqz

    # Sample rate and desired cutoff frequencies (in Hz).
    fs = 1.0/dt
    lowcut = fs/2000.0
    highcut = fs/5.1
    print dt,fs	
	 
    # Filter a noisy signal.
    plt.figure(3)
    plt.clf()
    plt.plot(time, h22_imag, label='With out applying filter')
#    plt.plot(time, h22_imag_SXS, label='SXS')
	
    h22_imag_filter = butter_bandpass_filter(h22_imag, lowcut, highcut, fs, order=3)
    plt.plot(time, h22_imag_filter, label='Filtered signal')
    plt.xlabel('time (seconds)')
   # plt.hlines([-a, a], 0, T, linestyles='--')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='lower left')

    plt.show()
