import matplotlib.pyplot as plt
import numpy as np
import plotsettings

def f(x):
	fx=np.exp(-x**2)
	return fx

t = np.linspace(-100, 100.0, 10000)
dt = t[2]-t[1]
ft=f(t)


F = np.fft.fft(ft)
freq = np.fft.fftfreq(len(ft),d=dt)

plt.plot(freq, abs(F))
plt.show()
