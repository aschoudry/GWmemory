import matplotlib.pyplot as plt
import numpy as np
import plottools
# >>>
fig,ax = plt.subplots()
x = np.linspace(0,1,100)
y = 1-x + 0.02*(2*np.random.random(len(x))-1)
ax.plot(x,y)
ax_zoom = plottools.zoom_axes(fig,ax,[0.1,0.2],[0.8,0.9],[0.6,0.9],[0.6,0.9])
ax_zoom.plot(x,y)
plt.show()
