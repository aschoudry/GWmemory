import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.integrate import odeint
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)

#memory expressions from Favata's paper
def H0mem(theta):
	H0 = (1.0/96.0)*(np.sin(theta)**2)*(17.0 + (np.cos(theta)**2))
	return H0

def H1mem(theta, eta):
	H1 = (np.sin(theta)**2)*(-(354241.0/2064384.0)-(62059.0/1032192.0)*(np.cos(theta)**2)-(4195.0/688128.0)*(np.cos(theta)**4) + \
		  ((15607.0/73728.0) + (9373.0/36864.0)*(np.cos(theta)**2) + (215.0/8192.0)*(np.cos(theta)**4))*eta)
	return H1

def H2mem(theta, eta):
	H2 = (np.sin(theta)**2)*(-(3968456539.0/9364045824.0) + (570408173.0/4682022912.0)*(np.cos(theta)**2) + \
		 (122166887.0/3121348608.0)*(np.cos(theta)**4) + (75601.0/15925248.0)*(np.cos(theta)**6) + (-(7169749.0/18579456.0)-\
		 (13220477.0/18579456.0)*(np.cos(theta)**2) - (1345405.0/6193152.0)*(np.cos(theta)**4) - (25115.0/884736.0)*(np.cos(theta)**6))*eta \
			+ ((10097.0/147456.0) + (5179.0/36864.0)*(np.cos(theta)**2) + (44765.0/147456.0)*(np.cos(theta)**4) + (3395.0/73728.0)*(np.cos(theta)**6))*eta**2)

	return H2

def H2p5mem(theta, eta):
	H2p5 = -5*(np.pi/21504.0)*(1.0- 4.0*eta)*(np.sin(theta)**2)*(509.0 + 472.0*(np.cos(theta)**2) + 39.0*(np.cos(theta)**4))
	return H2p5

def H3mem(theta, eta):
	H3 = (np.sin(theta)**2)*( -(69549016726181.0/46146017820672.0) + (6094001938489.0/23073008910336.0)*(np.cos(theta)**2) - (1416964616993.0/15382005940224.0)*(np.cos(theta)**4) \
		-(2455732667.0/78479622144.0)*(np.cos(theta)**6) - (9979199.0/2491416576.0)*(np.cos(theta)**8) + ((1355497856557.0/149824733184.0) - (3485.0*np.pi**2/9216) + \
		(-(3769402979.0/4682022912.0) - (205.0*np.pi**2/9216.0))*(np.cos(theta)**2) + (31566573919.0/49941577728.0)*(np.cos(theta)**4) + (788261497.0/3567255552.0)* \
		(np.cos(theta)**6) + (302431.0/9437184.0)*(np.cos(theta)**8))*eta + ((5319395.0/28311552.0) - (24019355.0/99090432.0)*(np.cos(theta)**2) - (4438085.0/3145728.0)* \
		(np.cos(theta)**4) - (3393935.0/7077888.0)*(np.cos(theta)**6) - (7835.0/98304.0)*(np.cos(theta)**8))*eta**2 + ((1433545.0/63700992.0) + (752315.0/15925248.0)* \
		(np.cos(theta)**2) + (129185.0/2359296.0)*(np.cos(theta)**4)+ (389095.0/1179648.0)*(np.cos(theta)**6)+ (9065.0/131072.0)*(np.cos(theta)**8))*eta**3)

	return H3


def dx_by_dt(x, M, eta, PN_order):
		
	gammaE = 0.57721
	C0 = 1.0
	C1 = -(743.0/336.0) -(11.0/4.0)*eta
	C1p5 = 4*np.pi
	C2 = (34103.0/18144.0) + (13661.0/2016.0)*eta + (59.0/18.0)*eta**2
	C2p5 = np.pi*( -(4159.0/672.0) - (189.0/8.0)*eta)
	C3 = (16447322263.0/139708800.0) + (16.0/3.0)*np.pi**2 - (856.0/105.0)*(2*gammaE + np.log(16.0*x)) + (-(56198689.0/217728.0) + (451.0/48.0)*np.pi**2)*eta + \
		(541.0/896.0)*eta**2 -(5605.0/2592.0)*(eta**3)
	C7p5 = np.pi*( -(4415.0/4032.0) + (358675.0/6048.0)*eta + (91495.0/1512.0)*eta**2)
	
	if PN_order==0:
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0)
	if PN_order==1:
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x)
	if PN_order==1.5:
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5))
	if PN_order==2:
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5) + C2*pow(x,2))
	if PN_order==2.5:
		 dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5) + C2*pow(x,2) + C2p5*pow(x,2.5)) 		 
	if PN_order==3:
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5) + C2*pow(x,2) + C2p5*pow(x,2.5) + C3*pow(x,3))
	if PN_order==3.5:	
		dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5) + C2*pow(x,2) + C2p5*pow(x,2.5) + C3*pow(x,3) + C7p5*pow(x,3.5)) 

	return dx_bydt




def Inegrate(x0, dt, nsteps, eta, M, PN_order):
	
	x1 = x0
	X=np.array([])

	k=0	
	for i in range(nsteps):
		kx1 = dx_by_dt(x1, M, eta, PN_order)*dt
		
		kx2 = dx_by_dt(x1 + 0.5*kx1 , M, eta, PN_order)*dt

		kx3 = dx_by_dt(x1 + 0.5*kx2, M, eta, PN_order)*dt

		kx4 = dx_by_dt(x1 + kx3, M, eta, PN_order)*dt

		x2 = x1 + ((kx1 + 2*kx2 + 2*kx3 + kx4)/6.0)

		x1 = x2

		X=np.append(X,x1)
		k+=1	

	return X

def h_plus_mem(theta, eta, M, R, x0, dt, nsteps, PN_order):

	x= Inegrate(x0, dt, nsteps, eta, M, PN_order)
	
	A = (2.0*eta*M*x/R)
	B0 = H0mem(theta)
	B1 = H1mem(theta, eta)
	B2 = H2mem(theta, eta)
	B2p5 = H2p5mem(theta, eta)
	B3 = H3mem(theta, eta)
	
	if PN_order==0:
		h_mem = A*(B0)
	if PN_order==1:
		h_mem = A*(B0 + B1*x) 
	if PN_order==2:
		h_mem = A*(B0 + B1*x + B2*pow(x,2))
	if PN_order==2.5:
		h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5))
	if PN_order==3:
		h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))
	
	return h_mem

def hx_plus_mem(theta, eta, M, R, x, dt, PN_order):

	A = (2.0*eta*M*x/R)
	B0 = H0mem(theta)
	B1 = H1mem(theta, eta)
	B2 = H2mem(theta, eta)
	B2p5 = H2p5mem(theta, eta)
	B3 = H3mem(theta, eta)
	
	if PN_order==0:
		h_mem = A*(B0)
	if PN_order==1:
		h_mem = A*(B0 + B1*x) 
	if PN_order==2:
		h_mem = A*(B0 + B1*x + B2*pow(x,2))
	if PN_order==2.5:
		h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5))
	if PN_order==3:
		h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))
	
	return h_mem

def phase(x, x0, eta, PN_order):

	gammaE = 0.57721

	A = -(1.0/32.0/eta/pow(x,5.0/2.0))
	B0 = 1
	B1 = (3715.0/1008.0 + 55.0*eta/12.0)
	B1p5 = -10*np.pi
	B2 = (15293365.0/1016064.0 + 27145.0*eta/1008.0 + 3085.0*eta**2/144.0)
	B2p5 = np.pi*np.log(x/x0)*(38645.0/1344.0 - 65.0/16.0*eta)
	B3 = (12348611926451.0/18776862720.0 - 160.0*np.pi**2/3 - (856.0/21.0)*(2*gammaE + np.log(16.0*x)) + (-15737765635.0/12192768.0 + 2255.0/48.0*np.pi**2)*eta \
		+ (76055.0/6912.0)*eta**2 - (127825.0/5184.0)*eta**3)


	if PN_order==0:
		ph = A*(B0)
	if PN_order==1:
		ph = A*(B0 + B1*x) 
	if PN_order==2:
		ph = A*(B0 + B1*x + B2*pow(x,2))
	if PN_order==2.5:
		ph = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5))
	if PN_order==3:
		ph = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))

	return ph
	 
t0 = -7000.0
time = np.linspace(t0, 0.0, 10000)
dt=time[1]-time[0]
M=1
R=1
theta = np.pi/2.0
eta=1.0/4.0
PN_order = 2.5
nsteps = len(time)
print len(time)
x0=pow(-5.00*M/(256.0*time[0]*eta), 1.0/4.0)

x1 = Inegrate(x0, dt, nsteps, eta, M, 1)
x2 = Inegrate(x0, dt, nsteps, eta, M, 2)
x3 = Inegrate(x0, dt, nsteps, eta, M, 2.5)

plt.plot(time, x1)
plt.plot(time, x2)
plt.plot(time, x3)
plt.show()


#
hp1PN = hx_plus_mem(theta, eta, M, R, x1, dt, 1)*np.exp(1j*phase(x1, x0, eta, 1))
hp2PN = hx_plus_mem(theta, eta, M, R, x2, dt, 2)*np.exp(1j*phase(x2, x0, eta, 2))
hp2p5PN = hx_plus_mem(theta, eta, M, R, x3, dt, 2.5)*np.exp(1j*phase(x3, x0, eta, 2.5))

plt.plot(time, abs(hp1PN))
#plt.plot(time, abs(hp2PN))
plt.plot(time, abs(hp2p5PN))

#plt.plot(time, hp)
plt.show()


#take FFT
hpf1 = np.fft.fft(hp1PN)
#hpf2 = np.fft.fft(hp2PN)
hpf3 = np.fft.fft(hp2p5PN)
freq = np.fft.fftfreq(len(hpf1),d=dt)

# Fourier tranform using SPA
#Fdot = (3.0/2.0)*(np.sqrt(x)/np.pi)*dx_by_dt(x, M, eta, PN_order)

#hpf_spa = abs(hp)/np.sqrt(Fdot)
#f = pow(x, 3.0/2.0)/(np.pi)



# Make frequency domain plot 
plt.loglog(freq, abs(hpf1*dt), 'r')
#plt.loglog(freq, abs(hpf2*dt))
plt.loglog(freq, abs(hpf3*dt), 'k')
plt.show()

#plt.loglog(f-0.0015, abs(hpf_spa))
plt.show()

########### Check dxdt growth rate ##############
#dxdt1 = dx_by_dt(x, M, eta, 1.5)
#dxdt2 = dx_by_dt(x, M, eta, 2.5)
#dxdt3 = dx_by_dt(x, M, eta, 3)

#plt.plot(time, dxdt1)
#plt.plot(time, dxdt2)
#plt.plot(time, dxdt3)

#plt.show()
#################################################

#plt.plot(x, Fdot)
#plt.show()

#shifted FFT
#take FFT

#hp = hp[8000:9500]
#hpfF = np.fft.fft(hp)
#freqF = np.fft.fftfreq(len(hpfF),d=dt)

#plt.loglog(freqF, abs(hpfF*dt))
#plt.show()




