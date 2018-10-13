import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

t = np.linspace(0, 1, num=20)
x=pow(t,3)

y = x
y_int = integrate.cumtrapz(y, x, initial=0)
plt.plot(x, y_int, 'ro', x,  0.5 * x**2, 'b-')
#plt.plot(x, 'ro')
plt.show()


