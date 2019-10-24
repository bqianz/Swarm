# from __future__ import division # true divide of arrays
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

a = 1
c = 2
n = 50

temp_mesh = np.linspace(0, 2*np.pi, n+1)
temp = temp_mesh[:-1]

def find_neighbours(current,n):

    i= current[0]
    j = current[1]

    neighbours = []

    i_prev = i-1
    i_next = i+1
    if i == 0:
        i_prev = n-1
    elif i == n-1:
        i_next = 0
    
    j_prev = j-1
    j_next = j+1
    if j == 0:
        j_prev = n-1
    elif j == n-1:
        j_next = 0

    return [(i_prev,j),(i_next,j),(i,j_prev),(i,j_next)]

def euclidean_dist(u,v,temp,c,a):
    x1, y1, z1 = tor2cart(temp[u[0]],temp[u[1]],c,a)
    x2, y2, z2 = tor2cart(temp[v[0]],temp[v[1]],c,a)

    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def Dijkstra(n,start,end,temp):

    vertices = [(a,b) for a in range(n) for b in range(n)]
    dist = {vertex: math.inf for vertex in vertices}
    prev = {vertex: None for vertex in vertices}

    dist[start] = 0
    current = start;

    while vertices and current != end:
        current = min(vertices, key=lambda vertex: dist[vertex])
        vertices.remove(current)
        print(len(vertices))
        for v in find_neighbours(current,n):
            if v in vertices:
                alt = dist[current] + euclidean_dist(current, v, temp, c, a)
                if alt < dist[v]:
                    dist[v] = alt
                    prev[v] = current

    current = end
    th_g = [temp[end[0]]]
    ph_g = [temp[end[1]]]
    while current != start:
        th_g.append(temp[current[0]])
        ph_g.append(temp[current[1]])
        current = prev[current]
    
    th_g.append(temp[start[0]])
    ph_g.append(temp[start[1]])
    th_g = np.asarray(th_g)
    ph_g = np.asarray(ph_g)
    return th_g, ph_g

def tor2cart(theta, phi, c, a):
    x = (c + a*np.cos(theta)) * np.cos(phi)
    y = (c + a*np.cos(theta)) * np.sin(phi)
    z = a * np.sin(theta)
    return x, y, z

def periodic_op_sc(a,b,op_func,period):
    if abs(b-a) > period / 2:
        if a < b:
            a = a + period
        else:
            b = b + period
    return op_func(a,b)

def periodic_operation(arr1, arr2, op_func, period):
    result = np.copy(arr1)
    for i in range(arr1.size):
        result[i] = periodic_op_sc(arr1[i],arr2[i],op_func,period)
    return result

def periodic_caliberate(arr, period):
    result = np.copy(arr)
    for i in range(arr.size):
        result[i] = float(Decimal(arr[i]) % Decimal(period))
    return result

def draw_with_vertices(x,y,z,xv,yv,zv):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, rstride=5, cstride=5, edgecolors='none', alpha = 0.3)
    ax.view_init(36, 26)
    
    # draw 
    ax.scatter(xv, yv, zv, c='r')
    return fig, ax

def draw_torus(x,y,z):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, rstride=1, cstride=1, edgecolors='none', alpha = 0.3)
    ax.view_init(36, 26)
    
    return fig, ax

def draw_colour(x,y,z,c_array):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, facecolors = c_array)
    ax.view_init(36, 26)
    
    return fig, ax

def general_segment(a,b,n):
    # return shortest path between a and b with mod n
    if a==b:
        return []
    half = int(n/2)
    s = np.sign(b-a)
    temp = np.arange(a,b,s)
    if temp.size > half:
        return list(np.arange(a,b-s*n,-s) % n)
    else:
        return list(temp % n)

def both_segments(a,b,n):
    if a==b:
        return [[],[]]
    half = int(n/2)
    s = np.sign(b-a)
    temp = np.arange(a,b,s)

    return [list(temp), list(np.arange(a,b-s*n,-s) % n)]

def initial_path(start,end):
    half = int(n/2)

    th_a = start[0]
    th_b = end[0]
    ph_a = start[1]
    ph_b = end[1]
    

    # two segment
    ph_seg= general_segment(ph_a,ph_b,n)
    th_seg1, th_seg2 = both_segments(th_a,th_b,n)

    if abs(th_a - half) < abs(th_b - half):
        # th_a closer to pi
        th1 = [th_a for i in ph_seg] + th_seg1
        ph1 = ph_seg + [ph_b for i in th_seg1]

        th2 = [th_a for i in ph_seg] + th_seg2
        ph2 = ph_seg + [ph_b for i in th_seg2]
    else:
        # th_b closer to pi
        th1 = th_seg1 + [th_b for i in ph_seg]
        ph1 = [ph_a for i in th_seg1] + ph_seg

        th2 = th_seg2 + [th_b for i in ph_seg]
        ph2 = [ph_a for i in th_seg2] + ph_seg

    # add end point
    th1 = th1 + [th_b]
    ph1 = ph1 + [ph_b]
    th2 = th2 + [th_b]
    ph2 = ph2 + [ph_b]

    # three segments: crossing theta = pi
    th_seg1 = general_segment(th_a, half, n)
    th_seg2 = general_segment(half, th_b, n)

    th3 = th_seg1 + [half for i in ph_seg] + th_seg2 + [th_b]
    ph3 = [ph_a for i in th_seg1] + ph_seg + [ph_b for i in th_seg2] + [ph_b]

    return temp[np.array(th1)], temp[np.array(ph1)], temp[np.array(th2)], temp[np.array(ph2)], temp[np.array(th3)], temp[np.array(ph3)]

def curve_length(u,v):
    # calculate length from a parametrization r(s) where s is from 0 to 1
    N = th.size - 1
    delta = 1/N
    delta2 = 2/N
    period = 2 * math.pi
    
    # reminder: fix periodic operations
    v_diff_0 = periodic_op_sc(th[1], th[0], operator.sub, period)
    u_diff_0 = periodic_op_sc(ph[1], ph[0], operator.sub, period)
    u_diff_end = periodic_op_sc(ph[-1], ph[-2], operator.sub, period)
    v_diff_end = periodic_op_sc(th[-1], th[-2], operator.sub, period)

    sum = math.sqrt( ( a * v_diff_0 / delta )**2 + ( (c + a*math.cos(v[0])) * u_diff_0 / delta )**2)
    sum += math.sqrt( ( a * v_diff_end / delta )**2 + ( (c + a*math.cos(v[0])) * u_diff_end / delta )**2)
    sum = sum/2

    for i in range(1,N):
        v_diff_i = periodic_op_sc(v[i+1], v[i-1], operator.sub, period)
        u_diff_i = periodic_op_sc(u[i+1], u[i-1], operator.sub, period)

        sum += math.sqrt( ( a * v_diff_i / delta2 )**2 + ( (c + a*math.cos(v[i])) * u_diff_i / delta2 )**2)

    return sum * delta

def ex1():
    var = 30
    n = 100
    c, a = 2, 1
    start = (0,0)
    end = (var,0)
    expected_length = 2 * math.pi * a * var/100

    return n, c, a, start, end, expected_length

def functional_iteration(th,ph):
    period = 2 * math.pi

    ph_diff = periodic_operation(ph[2:], ph[0:-2], operator.sub, period)
    th_diff = periodic_operation(th[2:], th[0:-2], operator.sub, period)
    ph_sum = periodic_operation(ph[2:], ph[0:-2], operator.add, period)
    th_sum = periodic_operation(th[2:], th[0:-2], operator.add, period)

    new_ph = np.copy(ph)
    frac1 = np.divide(a * np.sin(th), c + a * np.cos(th))
    new_ph[1:-1] = ph_sum/2 + np.multiply( frac1[1:-1], np.multiply( ph_diff, th_diff )) / 4

    new_th = np.copy(th)
    frac2 = np.multiply(np.sin(th)/a, c + a*np.cos(th))
    new_th[1:-1] = th_sum/2 + np.multiply( frac2[1:-1] , np.multiply(ph_diff,ph_diff) ) / 8

    return periodic_caliberate(new_th, period), periodic_caliberate(new_ph, period)

def iterate(th,ph):
    curvelength = curve_length(ph,th)
    count = 0
    norm = math.inf
    while norm > 10 ** (-3) and count < 250:
        th_new, ph_new = functional_iteration(th,ph,a,c)
        new_curve_length = curve_length(ph_new,th_new)
        norm = abs(curvelength - new_curve_length)/curvelength
        curvelength = new_curve_length
        th = th_new
        ph = ph_new
        count += 1
    return curvelength

def distance(start, end, n, temp):
    # (theta_index, ph_index)

    if start == end:
        return 0

    th1, ph1, th2, ph2, th3, ph3 = initial_path(n,start,end,temp)

    return min([iterate(th1,ph1), iterate(th2,ph2), iterate(th3,ph3)])

def distance(start,end):
    if start == end:
        return 0

    th1, ph1, th2, ph2, th3, ph3 = initial_path(n,start,end,temp)

def store_grid(n,temp):
    matfile = 'grid.mat'
    h = 2 * math.pi / n

    d = np.zeros([n, n, n], temp.dtype)
    d_th = np.zeros([n, n, n], temp.dtype)
    d_ph = np.zeros([n, n, n], temp.dtype)


    start = process_time()
    # generate d
    for i in range(n):
        start = (i,0)
        for j in range(n):
            for k in range(n):
                end = (j,k)
                d[i,j,k] = distance(start,end,n,temp)
            print("i,j = " + str(i) + ", " + str(j))
            print(process_time())
    
    # # generate d_th
    # for j in range(n):
    #     for i in range(1,n-1):
    #         d_th[i,j,:] = (d[i+1,j,:] - d[i-1,j,:]) / (2*h)
    #     d_th[0,j,:] = (d[1,j,:] - d[n-1,j,:]) / (2*h)
    #     d_th[n-1,j,:] = (d[0,j,:] - d[n-2,j,:]) / (2*h)

    # # generate d_ph
    # for i in range(n):
    #     for k in range(1,n-1):
    #         d_ph[i,:,k] = (d[i,:,k-1] - d[i,:,k+1]) / (2*h)
    #     d_ph[i,:,0] = (d[i,:,n-1] - d[i,:,1]) / (2*h)
    #     d_ph[i,:,n-1] = (d[i,:,n-2] - d[i,:,0]) / (2*h)

    scipy.io.savemat(matfile, mdict={'distance': d, 'd_theta': d_th, 'd_phi': d_ph})
    return matfile

def array_expand(instance, n):
    expanded = np.zeros([n+2,n+2])
    expanded[1:-1,1:-1] = instance

    expanded[-1,1:-1] = instance[0,:]
    expanded[0,1:-1] = instance[-1,:]
    expanded[1:-1, -1] = instance[:,0]
    expanded[1:-1,0] = instance[:,-1]

    expanded[0,0] = instance[-1,-1]
    expanded[0,-1] = instance[-1,0]
    expanded[-1,0] = instance[0,-1]
    expanded[-1,-1] = instance[0,0]

    return expanded

def recaliberate(array):
    period = 2 * math.pi
    for i in range(len(array)):
        temp = array[i]
        while temp < 0:
            temp = temp + period
        while temp >= period:
            temp = temp - period
        array[i] = temp
    return array


if __name__ == "__main__":

    # plotting
    temp_mesh = np.linspace(0, 2*np.pi, n+1)
    temp = temp_mesh[:-1]

    store_grid(n,temp)
