import numpy as np

def tor2cart(u, v, c, a):
    x = (c + a*np.cos(v)) * np.cos(u)
    y = (c + a*np.cos(v)) * np.sin(u)
    z = a * np.sin(v)
    return x, y, z




