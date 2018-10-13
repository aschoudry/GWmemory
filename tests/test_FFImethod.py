import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import fft

def f(x):
	return (np.cos(2*np.pi*x) + 12*np.cos(10*np.pi*x) + 10*np.cos(12*np.pi*x))

def compute_fourierTransform(t0, tf, dt, f):
	t=np.arange(t0, tf, dt)
	#Compute Fourier transform by numpy's FFT function
	g=np.fft.fft(f)
	#frequency normalization factor is 2*np.pi/dt
	w = np.fft.fftfreq(f.size)*2*np.pi/dt
	
	#In order to get a discretisation of the continuous Fourier transform
	#we need to multiply g by a phase factor
	#g*=dt*np.exp(-complex(0,1)*w*t0)/(np.sqrt(2*np.pi))

	return w, g 

def compute_integration_usingFFI(t0, tf, dt, w_cut, f):
	
	t=np.arange(t0, tf, dt)
	#Compute Fourier transform by numpy's FFT function
	g=np.fft.fft(f)
	#frequency normalization factor is 2*np.pi/dt
	w = np.fft.fftfreq(f.size, dt)*2*np.pi
	
	#In order to get a discretisation of the continuous Fourier transform
	#we need to multiply g by a phase factor
	g*=dt*np.exp(-complex(0,1)*w*t0)/(np.sqrt(2*np.pi))

	w[0]=w[1]
	dw=w[2]-w[1]

	for i in range(len(w)):
		if abs(w[i])<=w_cut:
			w[i]=np.sign(w[i])*w_cut

	invFT=np.fft.ifft(-1j*g/w)*dw/np.sqrt(2.0*np.pi)*len(t)

	return invFT

t0=-1.0
tf=1.0
dt=1.0/(12*2)

print 'fmax is', np.pi/dt

t=np.arange(t0,tf,dt)
y=compute_integration_usingFFI(t0, tf, dt, 0.0 ,f(t))
z=np.cumsum(f(t))*dt + y[0]
#plt.plot(t,f(t), 'k')
plt.plot(t,y, 'g')
plt.plot(t,z, 'r--')
plt.show()

	


