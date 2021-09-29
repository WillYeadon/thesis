import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np

radiusX = radiusY = np.arange(0, 2*np.pi, 0.01)  
radiusX, radiusY = np.meshgrid(radiusX, radiusY) 

x = radiusX*np.cos(radiusY)
y = radiusX*np.sin(radiusY)
z = (radiusX/2)**2

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
ax.plot_surface(x, y, z, cmap = cm.jet)

ax.view_init(elev=35., azim=315)

#ax.set_xticks(np.arange(-5, 6, step=2.5))
ax.set_xticklabels(np.arange(-7.5, 7.5, step = 2.5), fontsize = 25)
ax.set_yticklabels(np.arange(-7.5, 7.5, step = 2.5), fontsize = 25)
ax.set_zticklabels(np.arange(-10, 60, step=10), fontsize = 25)

ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)
ax.set_zlim(0, 13)

ax.set_xlabel('X', fontsize=30, labelpad=30)
ax.set_ylabel('Y', fontsize=30, labelpad=30)
ax.set_zlabel('Z', fontsize=30, labelpad=10)

fig.tight_layout(rect=[-1, -1, 1, 1])
fig.show()
plt.savefig('elliptic-paraboloid.png', bbox_inches='tight', dpi=300)