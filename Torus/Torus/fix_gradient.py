import scipy.io
import math
import numpy as np

n = 50
h = 2 * math.pi / n
temp_mesh = np.linspace(0, 2*np.pi, n+1)
temp = temp_mesh[:-1]


matdata = scipy.io.loadmat('grid.mat')
d =  matdata['distance']

d_th = np.zeros([n, n, n], temp.dtype)
d_ph = np.zeros([n, n, n], temp.dtype)

# generate d_th
for j in range(n):
    for i in range(1,n-1):
        d_th[i,j,:] = (d[i+1,j,:] - d[i-1,j,:]) / (2*h)
    d_th[0,j,:] = (d[1,j,:] - d[n-1,j,:]) / (2*h)
    d_th[n-1,j,:] = (d[0,j,:] - d[n-2,j,:]) / (2*h)

# generate d_ph
for i in range(n):
    for k in range(1,n-1):
        d_ph[i,:,k] = (d[i,:,k-1] - d[i,:,k+1]) / (2*h)
    d_ph[i,:,0] = (d[i,:,n-1] - d[i,:,1]) / (2*h)
    d_ph[i,:,n-1] = (d[i,:,n-2] - d[i,:,0]) / (2*h)

scipy.io.savemat('grid.mat', mdict={'distance': d, 'd_theta': d_th, 'd_phi': d_ph})
