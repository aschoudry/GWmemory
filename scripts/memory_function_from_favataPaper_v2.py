import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

#memory terms from arxiv.org/pdf/1003.3486.pdf equation 18
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

def dx_by_dt(x, t):
	
	#M=1
	#eta=0.25
	
	gammaE = 0.57721
	C0 = 1.0
	C1 = -(743.0/336.0) -(11.0/4.0)*eta
	C1p5 = 4*np.pi
	C2 = (34103.0/18144.0) + (13661.0/2016.0)*eta + (59.0/18.0)*eta**2
	C2p5 = np.pi*( -(4159.0/672.0) - (189.0/8.0)*eta)
	C3 = (16447322263.0/139708800.0) + (16.0/3.0)*np.pi**2 - (856.0/105.0)*(2*gammaE + np.log(16.0*x)) + (-(56198689.0/217728.0) + (451.0/48.0)*np.pi**2)*eta + (541.0/896.0)*eta**2 \
		-(5605.0/2592.0)*eta**3
	C7p5 = np.pi*( -(4415.0/4032.0) + (358675.0/6048.0)*eta + (91495.0/1512.0)*eta**2)

	dx_bydt = (64.0/5.0)*(eta/M)*pow(x, 5)*(C0 + C1*x + C1p5*pow(x, 1.5) + C2*pow(x,2) + C2p5*pow(x,2.5) + C3*pow(x,3) + C7p5*pow(x,7.5)) 

	return dx_bydt



def h_plus_mem_leadingOrder(theta, eta, M, R, t):
	
	x = pow(-5.0*M/(eta*256.0*t), 1.0/4.0) 
#	x = odeint(dx_by_dt, x0, t)
	
	A = (2.0*eta*M*x/R)
	B0 = H0mem(theta)
	B1 = H1mem(theta, eta)
	B2 = H2mem(theta, eta)
	B2p5 = H2p5mem(theta, eta)
	B3 = H3mem(theta, eta)
	
	h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))
	
	return h_mem

def h_plus_mem(theta, eta, M, R, t, x0):
	
#	x = pow(-2.0*M/(256.0*t), 1.0/4.0) 
	x = odeint(dx_by_dt, x0, t)
	
	A = (2.0*eta*M*x/R)
	B0 = H0mem(theta)
	B1 = H1mem(theta, eta)
	B2 = H2mem(theta, eta)
	B2p5 = H2p5mem(theta, eta)
	B3 = H3mem(theta, eta)
	
	h_mem = A*(B0 + B1*x + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))
	
	return h_mem


t = np.arange(-9000.0, -2000.0, 0.1)
M=1.0
q=2.0
eta=q/pow(1.0+q,2)
x0=pow(-5.0*M/(256.0*t[0]*eta), 1.0/4.0)

hp_mem_leading =  h_plus_mem_leadingOrder(np.pi/2, eta, 1.0, 1.0, t)
hp_mem = h_plus_mem(np.pi/2, eta, 1.0, 1.0, t, x0)

plt.plot(t, hp_mem)
plt.plot(t, hp_mem_leading,'r--')
plt.xlim(-9000, -1000)
plt.show()
 	
f=open("/home/ashok/gravitational_wave_memory_project/data/PostnewtonianMemory_v2.txt","w") 
for i in range(len(t)):
	f.write("%E %E\n" % (t[i], hp_mem[i])) 
f.close() 



