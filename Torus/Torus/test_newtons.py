# test if newton's iteration on an initialized closest path between two points generates the least distance path

import numpy as np
import scipy.io
from scipy import interpolate
import math
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.path import Path
from numpy import linalg as LA
import operator
from decimal import Decimal
from time import process_time

from newtonpath import NewtonPath
from Torus import initial_path


a = 1
c = 2

def it_newtons(u,v):
    count = 0
    Y = make_Y_intial(u,v)
    l = 4 * length(u)
    F = np.zeros(l)

    # continue here: specify stopping condition 
    while count < 50:
        J = make_Jacobian(Y)
        dY = LA.lstsq(J,F)
        sn = scaled_norm(dY)
        mu = determine_mu(sn)
        Y += mu * dY
        count += 1
    return curvelength

# find where is the code for drawing path!!
 

# generate two points on the Torus

# generate initial path between points on Torus

# plot Torus, two points, initial path

# newtons iteration on path, replot path at each iteration