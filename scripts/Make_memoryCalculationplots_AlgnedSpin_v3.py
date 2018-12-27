import numpy as np
import matplotlib.pyplot as plt
import plotsettings
from matplotlib.font_manager import FontProperties
from scipy.integrate import odeint

#memory expressions from Favata's paper(our results when spin included)
def H0mem(theta):
	H0 = (1.0/96.0)*(np.sin(theta)**2)*(17.0 + (np.cos(theta)**2))
	return H0

def H1mem(theta, eta):
	H1 = (np.sin(theta)**2)*(-(354241.0/2064384.0)-(62059.0/1032192.0)*(np.cos(theta)**2)-(4195.0/688128.0)*(np.cos(theta)**4) + \
		  ((15607.0/73728.0) + (9373.0/36864.0)*(np.cos(theta)**2) + (215.0/8192.0)*(np.cos(theta)**4))*eta)
	return H1
	
def H1p5(theta, eta, chiZa, chiZs):
	H1p5 = (np.sin(theta)**2)*(((3.0/1280.0) + (153.0/320.0)*np.sqrt(1.0 - 4*eta) - (3.0*eta/320.0) + (7.0/3840.0)*np.cos(theta)**2 + (9.0/320.0)*np.sqrt(1.0-4.0*eta)*np.cos(theta)**2\
		-(7.0*eta/960.0)*np.cos(theta)**2)*chiZa + ((153.0/320.0) + (3.0*np.sqrt(1.0 - 4.0*eta)/1280.0) - (3.0*eta/80.0) - (3.0*np.sqrt(1.0-4.0*eta)/320.0)*eta \
		+ (9.0/320.0)*np.cos(theta)**2  + (7.0*np.sqrt(1.0-4.0*eta)/3840.0)*np.cos(theta)**2 + (23.0/240.0)*eta*np.cos(theta)**2 -(7.0*np.sqrt(1.0-4.0*eta)/960.0)*eta\
		*np.cos(theta)**2)*chiZs) 	
	return H1p5

def H2mem(theta, eta):
	H2 = (np.sin(theta)**2)*((-(3968456539.0/9364045824.0) - (7169749.0*eta/18579456.0) + (10097.0/147456.0)*eta**2 + (570408173.0/4682022912.0)*np.cos(theta)**2\
		-(13220477.0/18579456.0)*eta*np.cos(theta)**2 + (5179.0/36864.0)*eta**2*np.cos(theta)**2 - (122166887.0/3121348608.0)*np.cos(theta)**4 - (1345405.0/6193152.0)*eta\
		*np.cos(theta)**4 + (44765.0/147456.0)*eta**2*np.cos(theta)**4 + (75601.0/15925248.0)*np.cos(theta)**6 - (25115.0*eta/884736.0)*np.cos(theta)**6 + (3395.0/73728.0)\
		*eta**2*np.cos(theta)**2) + (-(1385.0/4608.0))  )

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


def dx_by_dt(x, M, eta, PN_order, chiZa, chiZs):
		
	gammaE = 0.57721
	C0 = 1.0

	C1 = -(743.0/336.0) -(11.0/4.0)*eta

	C1p5 = 4*np.pi - (113.0/12.0)*np.sqrt(1.0- eta**2)*chiZa + (-(113.0/12.0) + (19.0/3.0)*eta)*chiZs

	C2 = (34103.0/18144.0) + (13661.0/2016.0)*eta + (59.0/18.0)*eta**2 + ((719.0/96.0)-30.0*eta)*chiZa**2 + ((-233.0/96.0)+ 10*eta)*chiZa**2 + \
		((719.0/48.0)*np.sqrt(1.0- eta**2)*chiZa*chiZs) + ((719.0/96.0) + (eta/24.0))*chiZs**2 - (233.0/48.0)*np.sqrt(1.0 - eta**2)*chiZa*chiZs + \
		((-233.0/96.0)-(7.0/24.0)*eta)*chiZs**2 

	C2p5 = np.pi*( -(4159.0/672.0) - (189.0/8.0)*eta) #+ np.sqrt(1.0-eta**2)*((-31319.0/1008.0)+(1159.0/24.0)*eta)*chiZa + \
		#((-31319.0/1008.0) + (22975.0/252.0)*eta - (79.0/3.0)*eta**2)*chiZs

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
	B0 = H0mem(theta)
	B1 = H1mem(theta, eta)
	B1p5 = H1p5(theta, eta, chiZa, chiZs)
	B2 = H2mem(theta, eta)
	B2p5 = H2p5mem(theta, eta)
	B3 = H3mem(theta, eta)
	
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
	
	h_mem[h_mem > 0.025]=0.25
	
	return h_mem

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SamMassRatio_differentSzSpinOnly/Memory_data/'
#import data
mass_ratio_vec = [0.43, 0.6, 0.8, 0.95]
filename_vec=['0p43', '0p600', '0p800','0p95']
i=0

for filename in filename_vec:
	
	datafile_hNRdot='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hdotNR.dat'
	datafile_hNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hNR.dat'
	datafile_hMemNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hMemNR.dat'

	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	timeNR, hdot_plus, hdot_cross = np.loadtxt(file_location+datafile_hNRdot, unpack=True)
	timeNR, h_plus, h_cross = np.loadtxt(file_location+datafile_hNR, unpack=True)


	#Making plots
	legend_size = 1
	fig = plt.figure()
	plt.suptitle('Memory Calculation for aligned Sz spin ='+str(mass_ratio_vec[i]) ,fontsize = 15)
	fontP = FontProperties()
	fontP.set_size('10.')


	legend = plt.legend(loc='best',prop={'size':legend_size})
	plt.subplot(2,2,1)
	plt.plot(timeNR, hdot_plus, 'r', label=r'$hdot_{+}$')
	plt.plot(timeNR, hdot_cross, 'k--', label=r'$hdot_{\times}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$hdot_{+,\times}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	plt.subplot(2,1,2)
	plt.plot(timeNR, h_plus, 'r', label=r'$h_{+}$')
	plt.plot(timeNR, h_cross, 'k--', label=r'$h_{\times}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$h_{+,\times}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	#legend = plt.legend(loc='best',prop={'size':legend_size})
	plt.subplot(2,2,2)
	plt.plot(timeNR, hmem, 'r', label=r'$h_{mem}$')
	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$h_{mem}$')
	plt.legend(loc=2)
	fontP.set_size('10.')

	plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_alignedSpin/'+filename+'.png')
	i+=1
	plt.close()


#data location
file_location ='/home/ashok/gravitational_wave_memory_project/data/SamMassRatio_differentSzSpinOnly/Memory_data/'


Spin_vec = [0.43, 0.60, 0.8, 0.95]
filename_vec=['0p43', '0p600', '0p800', '0p95']
tf_vec=[-60, -10, -10, -10.0]
ti_vec=[-8000, -8000, -8000, -8000]
idx_cut_vec = [80, 80, 80, 80, 80]


i=0
#Making plots
legend_size = 2
fig = plt.figure()
plt.title('Memory Calculation for Spinning case' ,fontsize = 15)
fontP = FontProperties()
fontP.set_size('20.')

legend = plt.legend(loc='best',prop={'size':legend_size})



for filename in filename_vec:
	
	datafile_hMemNR='rMPsi4_AlignedSpin_Sz_'+filename+'_q1dataClean_hMemNR.dat'
	timeNR, hmem, h_mem_plus = np.loadtxt(file_location+datafile_hMemNR, unpack=True)
	#Normalize to stich
	hmem*=17.0
	
	#generating Postnewtonian memory
	dt=timeNR[1]-timeNR[0]
	ti = ti_vec[i]
	tf=tf_vec[i]
	time_PN = np.arange(ti, tf, dt)
	nsteps=len(time_PN)+1
	
	q=1.0
	eta = q/pow(1.0+q,2)
	
	R=1.0
	M=1
	x0=pow(-5.0*M/(256.0*time_PN[0]*eta), 1.0/4.0)
	chiZa=0.0
	chiZs=Spin_vec[i]

 	hp_mem_PN = h_plus_mem(np.pi/2.0, eta, M, R, x0, dt, nsteps, 2, chiZa, chiZs)	
		

	time_PN=time_PN-tf
	dhp_mem_PN_dt=np.diff(hp_mem_PN)/dt
	dhp_mem_NR_dt=np.diff(hmem)/dt

	NRmax_idx = np.argmax(dhp_mem_NR_dt)
	dhp_mem_NR_dt=dhp_mem_NR_dt[:NRmax_idx]

	idx_cut = idx_cut_vec[i]
	dt=timeNR[1]-timeNR[0]
	Slope_hmem_NR=(hmem[idx_cut]-hmem[idx_cut-1])/dt
	idx=find_nearest_idx(dhp_mem_PN_dt, Slope_hmem_NR)
#	print idx, len(dhp_mem_PN_dt), q
	hmem_cut=hmem[idx_cut:]
	timeNR_cut=timeNR[idx_cut:]
	
#	print len(timeNR)
#	print timeNR_cut[0],time_PN[idx], timeNR_cut[0]-time_PN[idx]

	hp_mem_PN_cut=hp_mem_PN[:idx]
	time_PN_cut=np.linspace(ti-tf,timeNR_cut[0] , len(hp_mem_PN_cut))
	
	hmem_tot_cut=hmem_cut+hp_mem_PN_cut[-1]
	time_tot = np.append(time_PN_cut, timeNR_cut)
	hmem_tot = np.append(hp_mem_PN_cut, hmem_tot_cut)


	plt.title('Memory Calculation for Spinning Sz only = ' ,fontsize = 15)
	plt.plot(time_tot, hmem_tot, label=r'$hmem$ for $S_{z}$ = '+str(Spin_vec[i]))
	plt.plot(timeNR_cut, hmem_tot_cut, 'y--')

	plt.grid()
	plt.xlabel(r'$time$')
	plt.ylabel(r'$h_{mem}$')
	plt.legend(loc=2)
	fontP.set_size('20.')
	i+=1
	
	


plt.grid()
#plt.ylim(0,0.0007)
plt.xlabel(r'$time$')
plt.ylabel(r'$h_{mem}$')
plt.legend(loc=2)
fontP.set_size('10.')
plt.savefig('/home/ashok/gravitational_wave_memory_project/plots/MemoryPlot_alignedSpin/'+filename+'.pdf')
plt.show()
	
