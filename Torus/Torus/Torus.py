from __future__ import division # true divide of arrays
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.path import Path
from numpy import linalg as LA
import operator
from decimal import Decimal

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

def functional_iteration(th,ph,a,c):
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


def draw_base(x,y,z,xv,yv,zv):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, rstride=5, cstride=5, edgecolors='none', alpha = 0.3)
    ax.view_init(36, 26)
    
    # draw 
    ax.scatter(xv, yv, zv, c='r')
    return fig, ax

def general_segment(a,b,n):
    if a==b:
        return []
    half = int(n/2)
    s = np.sign(b-a)
    temp = np.arange(a,b,s)
    if temp.size > half:
        return list(np.arange(a,b-s*n,-s) % n)
    else:
        return list(temp % n)

def initial_path(n,start,end,temp):
    half = int(n/2)

    th_a = start[0]
    th_b = end[0]
    ph_a = start[1]
    ph_b = end[1]
    
    ph_seg= general_segment(ph_a,ph_b,n)
    th_seg = general_segment(th_a,th_b,n)

    if half in th_seg:
        # three segments - crossing theta = pi
        th_seg1 = general_segment(th_a, half, n)
        th_seg2 = general_segment(half, th_b, n)

        th = th_seg1 + [half for i in ph_seg] + th_seg2
        ph = [ph_a for i in th_seg1] + ph_seg + [ph_b for i in th_seg2]

    else:
        # two segment
        if abs(th_a - half) < abs(th_b - half):
            # th_a closer to pi
            th = [th_a for i in ph_seg] + th_seg
            ph = ph_seg + [ph_b for i in th_seg]
        else:
            # th_b closer to pi
            th = th_seg + [th_b for i in ph_seg]
            ph = [ph_a for i in th_seg] + ph_seg

    # add end point
    th = th + [th_b]
    ph = ph + [ph_b]

    th_values = temp[np.array(th)]
    ph_values = temp[np.array(ph)]

    return th_values, ph_values

def curve_length(th,ph):
    N = th.size - 1
    delta = 1/N
    delta2 = 2/N
    period = 2 * math.pi
    
    # reminder: fix periodic operations
    th_diff_0 = periodic_op_sc(th[1], th[0], operator.sub, period)
    ph_diff_0 = periodic_op_sc(ph[1], ph[0], operator.sub, period)
    ph_diff_end = periodic_op_sc(ph[-1], ph[-2], operator.sub, period)
    th_diff_end = periodic_op_sc(th[-1], th[-2], operator.sub, period)

    sum = math.sqrt( ( a * th_diff_0 / delta )**2 + ( (c + a*math.cos(th[0])) * ph_diff_0 / delta )**2)
    sum += math.sqrt( ( a * th_diff_end / delta )**2 + ( (c + a*math.cos(th[0])) * ph_diff_end / delta )**2)
    sum = sum/2

    for i in range(1,N):
        th_diff_i = periodic_op_sc(th[i+1], th[i-1], operator.sub, period)
        ph_diff_i = periodic_op_sc(ph[i+1], ph[i-1], operator.sub, period)

        sum += math.sqrt( ( a * th_diff_i / delta2 )**2 + ( (c + a*math.cos(th[i])) * ph_diff_i / delta2 )**2)

    return sum * delta

def ex1():
    var = 30
    n = 100
    c, a = 2, 1
    start = (0,0)
    end = (var,0)
    expected_length = 2 * math.pi * a * var/100

    return n, c, a, start, end, expected_length



if __name__ == "__main__":
    n = 100
    c, a = 2, 1

    # (theta_index, ph_index)
    start = (0,30)
    end = (0,10)

    expected_length = 0

    n, c, a, start, end, expected_length = ex1()


    # plotting
    temp_mesh = np.linspace(0, 2*np.pi, n+1)
    temp = temp_mesh[:-1]
    phi, theta = np.meshgrid(temp_mesh,temp_mesh)
    xx,yy,zz = tor2cart(theta,phi,c,a)
    xv,yv,zv = tor2cart(temp[np.array([start[0],end[0]])],temp[np.array([start[1],end[1]])],c,a)
    fig1, ax1 = draw_base(xx,yy,zz,xv,yv,zv)

    # get initial path
    th, ph = initial_path(n,start,end,temp)
    # th, ph = Dijkstra(n,start,end,temp)

    initial_curve_length = curve_length(th,ph)
    x,y,z = tor2cart(th,ph,c,a)
    ax1.plot(x,y,z,c='g')
    
    # smooth path with fixed point iteration
    count = 0
    norm = math.inf
    while norm > 0.0001 and count < 2000:
        th_new, ph_new = functional_iteration(th,ph,a,c)
        diff_th = LA.norm(th_new - th, np.inf)
        diff_ph = LA.norm(ph_new - ph, np.inf)
        norm= LA.norm([diff_th,diff_ph])
        th = th_new
        ph = ph_new
        count += 1

    fig2, ax2 = draw_base(xx,yy,zz,xv,yv,zv)
    x,y,z = tor2cart(th,ph,c,a)
    ax2.plot(x,y,z,c='g')

    print("convergence norm:")
    print(norm)

    print("initial curve length:")
    print(initial_curve_length)

    print("iterated curve length:")
    print(curve_length(th,ph))

    print("expected curve length:")
    print(expected_length)

    plt.show()