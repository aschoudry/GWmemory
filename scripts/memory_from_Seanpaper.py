import numpy as np
import matplotlib.pyplot as plt

def Omega(Omega0, k, t, tp, t0, tau):
	Omg=(pow(Omega0,4) + k*(np.arctan2(t-tp, tau)- np.arctan2(t0-tp, tau)))**(1.0/4.0)
	return Omg

Omega0=0
k=1
tp=0
t0=-1000
tau =1.0

time=np.arange(-1000, -10, 0.5)
hmem=np.array([])
for i in range(len(time)):
	hmem=np.append(hmem,Omega(Omega0, k, i, tp, t0, tau))
	

plt.plot(time,-hmem)
plt.show()
