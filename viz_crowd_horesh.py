

#
#
# the first command line argument is the name of the output file which
# from crowd_horesh
#
#
#



import sys
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

file = open(sys.argv[1], 'r')
Lines = file.readlines()

y=0
start = False
for line in Lines:
    line = line.strip()
    x = line.split(' ')[0][0:-1]
    if start:
        yprev = y
        y = int(line.split(' ')[1][0:-1])
        if y == 0:
            ydim = yprev+1
    if x == '0':
        start = True

xdim = int(Lines[-1].split(' ')[0][0:-1])+1

print(xdim,ydim)

xdata=np.zeros((xdim,ydim))
ydata=np.zeros((xdim,ydim))
zdata=np.zeros((xdim,ydim))

xcount = 0
ycount = 0
start = False
# Strips the newline character
for line in Lines:
    line = line.strip()
    x = line.split(' ')[0][0:-1]
    if x == '0':
        start = True
    if start:
        x = int(x)
        y = int(line.split(' ')[1][0:-1])
        z = float(line.split(' ')[2])
        xdata[x][y]= x
        ydata[x][y]= y
        zdata[x][y] = z


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.


# Plot the surface.
surf = ax.plot_surface(xdata, ydata, zdata, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
