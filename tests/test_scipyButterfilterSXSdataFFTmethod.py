from scipy.signal import butter, lfilter
import numpy as np
from scipy import integrate
from scipy import signal
from scipy import interpolate
import matplotlib.pyplot as plt
import h5py
from numpy import ceil, log2
from collections import deque

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

time_intrp=np.linspace(time[0], time[-1], len(time))
d2h22_real_intrp=interpolate.interp1d(time_intrp, d2h22_real_by_dt2)
d2h22_imag_intrp=interpolate.interp1d(time_intrp, d2h22_imag_by_dt2)

time=time_intrp
h22_real=h22_real_intrp(time)
h22_imag=h22_imag_intrp(time)
psi_tot=d2h22_real_by_dt2 + 1j*d2h22_imag_by_dt2
phase=-np.unwrap(np.angle(psi_tot))
instfreq=np.gradient(phase, time)



dt=time[2]-time[1]


#integrating the data second time with respect to time
def nextpow2(x):
    """returns the smallest power of two that is greater than or equal to the
    absolute value of x.

    This function is useful for optimizing FFT operations, which are
    most efficient when sequence length is an exact power of two.

    :Example:

    .. doctest::

        >>> from spectrum import nextpow2
        >>> x = [255, 256, 257]
        >>> nextpow2(x)
        array([8, 8, 9])

    """
    res = ceil(log2(x))
    return res.astype('int')  #we want integer values only but ceil gives float

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

fftsize=2**(nextpow2(len(time)))
hfreq = ((time[2]-time[1])/(fftsize/2.0 +1))*np.array(range(fftsize+1))
hfreq[fftsize/2+2:fftsize+1] = -hfreq[fftsize/2:1:-1]
hfreq[1] = 1e-30
if min(instfreq) > 0:
	fo = dsearchn(hfreq,min(instfreq)/5)
	fe = dsearchn(hfreq,-min(instfreq)/5)
else:
	fo = 1
        fe = fftsize

    psi_fft = (timeNR(2)-timeNR(1))*fft(psi_tot,fftsize)
    hdotfft = psi_fft/(timeNR(2)-timeNR(1));
    hdotfft(1:fo) = -hdotfft(1:fo)./hfreq(fo)/2/pi;
    hdotfft(fe:end) = -hdotfft(fe:end)./hfreq(fo)/2/pi;
    hdotfft(fo+1:fe-1) = -hdotfft(fo+1:fe-1)./hfreq(fo+1:fe-1)/2/pi;
    hdot_tot = ifft(hdotfft,fftsize);
    hdot_tot = hdot_tot(1:length(timeNR));
    hdot_plus = real(hdot_tot);
    hdot_cross = imag(hdot_tot);
    hfft = psi_fft/(timeNR(2)-timeNR(1));
    hfft(1:fo) = -hfft(1:fo)./hfreq(fo)./hfreq(fo)/4/pi/pi;
    hfft(fe:end) = -hfft(fe:end)./hfreq(fo)./hfreq(fo)/4/pi/pi;
    hfft(fo+1:fe-1) = -hfft(fo+1:fe-1)./hfreq(fo+1:fe-1)./hfreq(fo+1:fe-1)/4/pi/pi;
    h_tot = ifft(hfft,fftsize);
    h_tot = h_tot(1:length(timeNR));
    h_plus = real(h_tot);
    h_cross = imag(h_tot)
