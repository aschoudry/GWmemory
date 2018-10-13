import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

t = np.linspace(0, 1, num=20)
x=pow(t,3)

y = x*x
yp = np.gradient(y, x)
y_int = integrate.cumtrapz(yp, x, initial=0)
plt.plot(x, y_int, 'ro', x,  x**2, 'b-')
#plt.plot(x, 'ro')
plt.show()

