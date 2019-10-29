
import numpy as np

matfile = 'torus_data.mat'

c, a = 2, 1

temp = np.linspace(0, 2*np.pi, 31)
u,v = np.meshgrid(temp[:-1],temp[:-1])

def tor2cart(u,v):
    x = (c + a*np.cos(v)) * np.cos(u)
    y = (c + a*np.cos(v)) * np.sin(u)
    z = a * np.sin(v)
    return x, y, z

x,y,z = tor2cart(u,v)

scipy.io.savemat(matfile, mdict={'xyz': (x,y,z)})