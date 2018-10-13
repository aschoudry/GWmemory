import numpy as np
import matplotlib.pyplot as plt

def A(t, t0, t1, A1, A2, sigma0, sigma2):
	amp=(A1/2.0)*(1.0 + np.tanh((t-t0)/sigma0))*(1.0 + A2*np.exp(t/sigma2))*(1.0 + np.tanh(-(t-t1)/sigma1))
	return amp

def Phi(t, t1, tf, omega_i, omega_f, sigmaPhi):
	phi=omega_i*(t-t1) + ((omega_f - omega_i)/2.0)*(1.0 + sigmaPhi*np.log(np.cosh((t-tf)/sigmaPhi)))
	return phi	

def Omega(t):
	omega = omega_i + ((omega_f - omega_i)/2.0)*(1.0 + sigma)
	return omega

def h(t, t0, t1, A1, A2, sigma0, sigma2, tf, omega_i, omega_f, sigmaPhi):
	h= A(t, t0, t1, A1, A2, sigma0, sigma2)*np.exp(1j*Phi(t, t1, tf, omega_i, omega_f, sigmaPhi))
	return h.real

t0=-430.0
t1=0.0
tf=100
omega_i=0.2
omega_f=1.0
A1=0.002
A2=5.0
sigma0=10.0
sigma1=16.0
sigma2=80.0
sigmaPhi=80.0

print Phi(10, t1, tf, omega_i, omega_f, sigmaPhi)


t=np.arange(-500,100,0.5)
#plt.plot(t, Phi(t, t1, tf, omega_i, omega_f, sigmaPhi))
plt.plot(t,  h(t, t0, t1, A1, A2, sigma0, sigma2, tf, omega_i, omega_f, sigmaPhi))
plt.show()
