import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.integrate import odeint

#memory expressions from Favata's paper(our results when spin included)
def H0mem(theta, eta, chiZa, chiZs):
	H0 = (0.17708333333333334 + np.cos(theta)**2/96.)*np.sin(theta)**2
	return H0

def H1mem(theta, eta, chiZa, chiZs):
	H1 = (-0.17159646654885913 - (62059*np.cos(theta)**2)/1.032192e6 - (4195*np.cos(theta)**4)/688128. +\
		eta*(0.21168348524305555 + (9373*np.cos(theta)**2)/36864. + (215*np.cos(theta)**4)/8192.))*np.sin(theta)**2
	return H1
	
def H1p5(theta, eta, chiZa, chiZs):
	H1p5 = (chiZa*(0.00234375 + (153*np.sqrt(1 - 4*eta))/320. + (7*np.cos(theta)**2)/3840. + (9*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/320. + \
		eta*(-0.009375 - (7*np.cos(theta)**2)/960.)) +\
		chiZs*(0.478125 + (3*np.sqrt(1 - 4*eta))/1280. + (9*np.cos(theta)**2)/320. + (7*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/3840. + \
		eta*(-0.0375 - (3*np.sqrt(1 - 4*eta))/320. + (23*np.cos(theta)**2)/240. - (7*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/960.)))*np.sin(theta)**2
	return H1p5

def H2mem(theta, eta, chiZa, chiZs):
	H2 = ((-0.4237972147497257 + (570408173*np.cos(theta)**2)/4.682022912e9 + (122166887*np.cos(theta)**4)/3.121348608e9 +\
		75601*np.cos(theta)**6)/1.5925248e7 + eta*\
		(-0.38589660536885473 - (13220477*np.cos(theta)**2)/1.8579456e7 - (1345405*np.cos(theta)**4)/6.193152e6 -\
		(25115*np.cos(theta)**6)/884736.) + eta**2*\
		(0.06847466362847222 + (5179*np.cos(theta)**2)/36864. + (44765*np.cos(theta)**4)/147456. + \
		(3395*np.cos(theta)**6)/73728.) + chiZs**2*\
		(-0.3005642361111111 - (85*np.cos(theta)**2)/4608. + eta**2*(-0.027777777777777776 - np.cos(theta)**2/72.) + \
		eta*(0.028645833333333332 + np.cos(theta)**2/128.)) + \
		chiZa**2*(-0.3005642361111111 - (85*np.cos(theta)**2)/4608. + eta*(1.1875 + (7*np.cos(theta)**2)/96.)) +\
		chiZa*chiZs*((-1385*np.sqrt(1 - 4*eta))/2304. - (85*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/2304. + \
		eta*(np.sqrt(1 - 4*eta)/72. + (np.sqrt(1 - 4*eta)*np.cos(theta)**2)/144.)))*np.sin(theta)**2

	return H2

def H2p5mem(theta, eta, chiZa, chiZs):
	H2p5 =((-2545*np.pi)/21504. - (295*np.pi*np.cos(theta)**2)/2688. - (65*np.pi*np.cos(theta)**4)/7168. +\
		eta*((2545*np.pi)/5376. + (295*np.pi*np.cos(theta)**2)/672. + (65*np.pi*np.cos(theta)**4)/1792.) +\
		chiZs*(1.6932645006244684 + (11867*np.sqrt(1 - 4*eta))/3.612672e6 - (484979*np.cos(theta)**2)/2.408448e6 +\
		(163*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/112896. - (70645*np.cos(theta)**4)/2.064384e6 +\
		(55*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/1.204224e6 + \
		eta*(-2.6840451241260395 - (12343*np.sqrt(1 - 4*eta))/1.806336e6 + (2001211*np.cos(theta)**2)/1.354752e6 -\
		(1517*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/903168. + (24725*np.cos(theta)**4)/129024. - \
		(55*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/200704.) +\
		eta**2*(0.7254866573365457 - (3797*np.sqrt(1 - 4*eta))/150528. - (1072289*np.cos(theta)**2)/1.354752e6 - \
		(411*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/25088. - (16785*np.cos(theta)**4)/100352. + \
		(55*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/150528.)) + \
		chiZa*(0.0032848263003117913 + (8156279*np.sqrt(1 - 4*eta))/4.816896e6 + (163*np.cos(theta)**2)/112896. -\
		(484979*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/2.408448e6 + (55*np.cos(theta)**4)/1.204224e6 - \
		(70645*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/2.064384e6 + \
		eta**2*(-0.02522454294217687 - (411*np.cos(theta)**2)/25088. + (55*np.cos(theta)**4)/150528.) + \
		eta*(-0.006833169465702948 - (4105391*np.sqrt(1 - 4*eta))/1.0838016e7 - (1517*np.cos(theta)**2)/903168. + \
		(7446571*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/5.419008e6 - (55*np.cos(theta)**4)/200704. + \
		(58395*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/401408.)))*np.sin(theta)**2

	return H2p5

def H3mem(theta, eta, chiZa, chiZs):
	H3 =(-1.5071509423140987 + (6093995955001*np.cos(theta)**2)/2.3073008910336e13 - \
		(1416977186081*np.cos(theta)**4)/1.5382005940224e13 - (2455652411*np.cos(theta)**6)/7.8479622144e10 - \
		(9979199*np.cos(theta)**8)/2.491416576e9 +\
		eta**2*(0.18788881486896064 - (792654827*np.cos(theta)**2)/3.269984256e9 - \
		(9226880251*np.cos(theta)**4)/6.539968512e9 - (71269747*np.cos(theta)**6)/1.48635648e8 -\
		 (7835*np.cos(theta)**8)/98304.) + eta*\
		(9.047223133581763 - (3485*np.pi**2)/9216. - (41463325921*np.cos(theta)**2)/5.1502252032e10 -\
		(205*np.pi**2*np.cos(theta)**2)/9216. + (347235904277*np.cos(theta)**4)/5.49357355008e11 +\
		(788232313*np.cos(theta)**6)/3.567255552e9 + (302431*np.cos(theta)**8)/9.437184e6) +\
		eta**3*(0.02250363474940637 + (57932071*np.cos(theta)**2)/1.226244096e9 + \
		(29847079*np.cos(theta)**4)/5.44997376e8 + (8170691*np.cos(theta)**6)/2.4772608e7 + \
		(9065*np.cos(theta)**8)/131072.) + chiZs*\
		((4675*np.pi)/1152. - (9*np.sqrt(1 - 4*eta)*np.pi)/2048. + (275*np.pi*np.cos(theta)**2)/1152. - \
		(7*np.sqrt(1 - 4*eta)*np.pi*np.cos(theta)**2)/2048. + \
		eta*((-947*np.pi)/288. + (9*np.sqrt(1 - 4*eta)*np.pi)/512. - (91*np.pi*np.cos(theta)**2)/288. +\
		(7*np.sqrt(1 - 4*eta)*np.pi*np.cos(theta)**2)/512.)) + \
		chiZs**2*(-0.21864656418088882 + (113*np.sqrt(1 - 4*eta))/8192. + (3716075*np.cos(theta)**2)/3.3030144e7 + \
		(791*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/73728. + (113265*np.cos(theta)**4)/7.340032e6 + \
		 eta**3*(-0.057291666666666664 - (11*np.cos(theta)**2)/384.) + \
		eta*(5.739628746396019 - (33*np.sqrt(1 - 4*eta))/512. + (921469*np.cos(theta)**2)/4.128768e6 - \
		(77*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/1536. - (184955*np.cos(theta)**4)/2.752512e6) + \
		eta**2*(-2.302741883293031 + (19*np.sqrt(1 - 4*eta))/512. - (676133*np.cos(theta)**2)/2.064384e6 + \
		(133*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/4608. + (11745*np.cos(theta)**4)/458752.)) + \
		chiZa**2*(-0.21864656418088882 + (113*np.sqrt(1 - 4*eta))/8192. + (3716075*np.cos(theta)**2)/3.3030144e7 + \
		(791*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/73728. + (113265*np.cos(theta)**4)/7.340032e6 + \
		eta*(1.9677733769492498 - (113*np.sqrt(1 - 4*eta))/2048. - (8200699*np.cos(theta)**2)/8.257536e6 -\
		(791*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/18432. - (701315*np.cos(theta)**4)/5.505024e6) + \
		eta**2*(-3.8601616753472223 + (40501*np.cos(theta)**2)/18432. + (1075*np.cos(theta)**4)/4096.)) + \
		chiZa*((-9*np.pi)/2048. + (4675*np.sqrt(1 - 4*eta)*np.pi)/1152. - (7*np.pi*np.cos(theta)**2)/2048. + \
		(275*np.sqrt(1 - 4*eta)*np.pi*np.cos(theta)**2)/1152. + eta*((9*np.pi)/512. + (7*np.pi*np.cos(theta)**2)/512.) +\
		chiZs*(0.027587890625 - (14443855*np.sqrt(1 - 4*eta))/3.3030144e7 + (791*np.cos(theta)**2)/36864. + \
		(3716075*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/1.6515072e7 + \
		(113265*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/3.670016e6 + \
		eta**2*(0.2578125 + (11*np.sqrt(1 - 4*eta))/384. + (77*np.cos(theta)**2)/384. + \
		(11*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/768.) + \
		eta*(-0.1748046875 + (56422223*np.sqrt(1 - 4*eta))/8.257536e6 - (1253*np.cos(theta)**2)/9216. - \
		(440281*np.sqrt(1 - 4*eta)*np.cos(theta)**2)/1.376256e6 - (17415*np.sqrt(1 - 4*eta)*np.cos(theta)**4)/131072.))))*np.sin(theta)**2
	
	return H3


def dx_by_dt(x, M, eta, PN_order, chiZa, chiZs):
		
	gammaE = 0.57721
	C0 = 1.0

	C1 = -(743.0/336.0) -(11.0/4.0)*eta

	C1p5 = 4*np.pi - (113.0/12.0)*np.sqrt(1.0- eta)*chiZa + (-(113.0/12.0) + (19.0/3.0)*eta)*chiZs

	C2 = (34103.0/18144.0) + (13661.0/2016.0)*eta + (59.0/18.0)*eta**2 + ((719.0/96.0)-30.0*eta)*chiZa**2 + ((-233.0/96.0)+ 10*eta)*chiZa**2 + \
		((81.0/8.0)*np.sqrt(1.0- eta**2)*chiZa*chiZs) + ((719.0/96.0) + (eta/24.0))*chiZs**2 - (233.0/48.0)*np.sqrt(1.0 - eta**2)*chiZa*chiZs + \
		((-233.0/96.0)-(7.0/24.0)*eta)*chiZs**2 

	C2p5 = np.pi*( -(4159.0/672.0) - (189.0/8.0)*eta) + np.sqrt(1.0-eta**2)*((-31319.0/1008.0)+(1159.0/24.0)*eta)*chiZa + \
		((-31319.0/1008.0) + (22975.0/252.0)*eta - (79.0/3.0)*eta**2)*chiZs

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




def Inegrate(x0, dt, nsteps, eta, M, PN_order, chiZa, chiZs):
	
	x1 = x0
	X=np.array([])

	k=0	
	for i in range(nsteps):
		kx1 = dx_by_dt(x1, M, eta, PN_order, chiZa, chiZs)*dt
		
		kx2 = dx_by_dt(x1 + 0.5*kx1 , M, eta, PN_order, chiZa, chiZs)*dt

		kx3 = dx_by_dt(x1 + 0.5*kx2, M, eta, PN_order, chiZa, chiZs)*dt

		kx4 = dx_by_dt(x1 + kx3, M, eta, PN_order, chiZa, chiZs)*dt

		x2 = x1 + ((kx1 + 2*kx2 + 2*kx3 + kx4)/6.0)

		x1 = x2

		X=np.append(X,x1)
		k+=1		

	return X

def h_plus_mem(theta, eta, M, R, x0, dt, nsteps, PN_order, chiZa, chiZs):

	x= Inegrate(x0, dt, nsteps, eta, M, PN_order, chiZa, chiZs)
	
	A = (2.0*eta*M*x/R)
	B0 = H0mem(theta, eta, chiZa, chiZs)
	B1 = H1mem(theta, eta, chiZa, chiZs)
	B1p5 = H1p5(theta, eta, chiZa, chiZs)
	B2 = H2mem(theta, eta, chiZa, chiZs)
	B2p5 = H2p5mem(theta, eta, chiZa, chiZs)
	B3 = H3mem(theta, eta, chiZa, chiZs)
	
	if PN_order==0:
		h_mem = A*(B0)
	if PN_order==1:
		h_mem = A*(B0 + B1*x) 
	if PN_order==1.5:
		h_mem = A*(B0 + B1*x + B1p5*pow(x,1.5))	
	if PN_order==2:
		h_mem = A*(B0 + B1*x + B1p5*pow(x,1.5)+ B2*pow(x,2))
	if PN_order==2.5:
		h_mem = A*(B0 + B1*x + B1p5*pow(x,1.5) + B2*pow(x,2) + B2p5*pow(x,2.5))
	if PN_order==3:
		h_mem = A*(B0 + B1*x + B1p5*pow(x,1.5) + B2*pow(x,2) + B2p5*pow(x,2.5) + B3*pow(x,3))
	
	return h_mem

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx



