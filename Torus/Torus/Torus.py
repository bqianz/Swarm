from __future__ import division # true divide of arrays
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.path import Path
from numpy import linalg as LA

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


def functional_iteration(th,ph,a,c):
    ph_diff = ph[2:] - ph[0:-2]
    th_diff = th[2:] - th[0:-2]

    new_ph = np.copy(ph)

    term1 = np.divide(a * np.sin(th), c + a * np.cos(th))
    new_ph[1:-1] = (ph[2:] + ph[0:-2])/2 + np.multiply( term1[1:-1], np.multiply( ph_diff, th_diff )) / 4

    new_th = np.copy(th)
    term2 = np.multiply(np.sin(th)/a, c + a*np.cos(th))
    new_th[1:-1] = (th[2:] + th[0:-2])/2 + np.multiply( term2[1:-1] , np.multiply(ph_diff,ph_diff) ) / 8

    return new_th, new_ph


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

def initial_path(n,start,end,temp):
    half_th = int(n/2)

    th_s = start[0]
    th_e = end[0]
    ph_s = start[1]
    ph_e = end[1]
    
    if abs(th_s - th_e) < half_th:
        # go through half_th
        if th_s < th_e:
            th_seg1 = range(start[0],half_th)
            th_seg2 = range(half_th,end[0])
    else:
        # go through th = 0
        ph_seg = 0

    seg1 = range(start[0],half_th)
    seg2 = range(start[1],end[1])
    seg3 = range(half_th,end[0])

    th_ind = list(seg1) + [half_th for i in seg2] + list(seg3) + [end[0]]
    ph_ind = [start[1] for i in seg1] + list(seg2) + [end[1] for i in seg3] + [end[1]]
    th = temp[np.array(th_ind)]
    ph = temp[np.array(ph_ind)]
    return th, ph, th.size

def curve_length(th,ph):
    N = th.size - 1
    delta = 1/N
    delta2 = 2/N

    sum = math.sqrt( ( a * (th[1] - th[0]) / delta )**2 + ( (c + a*math.cos(th[0])) * (ph[1] - ph[0]) / delta )**2)
    sum += math.sqrt( ( a * (th[-1] - th[-2]) / delta )**2 + ( (c + a*math.cos(th[0])) * (ph[-1] - ph[-2]) / delta )**2)
    sum = sum/2

    for i in range(1,N):
        sum += math.sqrt( ( a * (th[i+1] - th[i-1]) / delta2 )**2 + ( (c + a*math.cos(th[i])) * (ph[i+1] - ph[i-1]) / delta2 )**2)

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
    start = (25,25)
    end = (0,0)

    expected_length = 0

    # n, c, a, start, end, expected_length = ex1()


    # plotting
    temp_mesh = np.linspace(0, 2*np.pi, n+1)
    temp = temp_mesh[:-1]
    phi, theta = np.meshgrid(temp_mesh,temp_mesh)
    xx,yy,zz = tor2cart(theta,phi,c,a)
    xv,yv,zv = tor2cart(temp[np.array([start[0],end[0]])],temp[np.array([start[1],end[1]])],c,a)
    fig1, ax1 = draw_base(xx,yy,zz,xv,yv,zv)

    # get initial path
    # th, ph, N = initial_path(n,start,end,temp)

    th, ph = Dijkstra(n,start,end,temp)

    print(th)
    print(ph)

    x,y,z = tor2cart(th,ph,c,a)
    ax1.plot(x,y,z,c='g')
    
    # smooth path with fixed point iteration
    count = 0
    norm = math.inf
    while norm > 0.0001:
        th_new, ph_new = functional_iteration(th,ph,a,c)
        diff_th = LA.norm(th_new - th, np.inf)
        diff_ph = LA.norm(ph_new - ph, np.inf)
        norm = LA.norm([diff_th,diff_ph])
        th = th_new
        ph = ph_new
        count += 1
        print(count)
    
    print("norm:")
    print(norm)


    fig2, ax2 = draw_base(xx,yy,zz,xv,yv,zv)
    x,y,z = tor2cart(th,ph,c,a)
    ax2.plot(x,y,z,c='g')

    print(th)
    print(ph)
    
    print("expected length:")
    print(expected_length)

    print("curve length:")
    print(curve_length(th,ph))

    plt.show()