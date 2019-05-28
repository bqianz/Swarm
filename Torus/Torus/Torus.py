from __future__ import division # true divide of arrays
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.path import Path

def find_neighbours(current,visited,n,temp):

    i= current[0]
    j = current[1]

    neighbours = []
    neighbours_th = []
    neighbours_ph = []

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

    for u in [[i_prev,j],[i_next,j],[i,j_prev],[i,j_next]]:
        if visited[u[0],u[1]] == 0:
            neighbours.append(u)
            neighbours_th.append(temp[u[0]])
            neighbours_ph.append(temp[u[1]])
    return neighbours_th, neighbours_ph, neighbours


def Dijkstra(n,start,end,temp):

    dist = np.full((n,n), n**2)
    prev = np.full((n,n,2),-1)
    visited = np.zeros((n,n))

    dist[start[0],start[1]] = 0
    visited[start[0],start[1]] = 1
    current = start;

    while np.prod(visited) == 0 and not np.array_equal(current,end):
        current_dist = dist[current[0],current[1]]
        _, _, neighbours = find_neighbours(current,visited,n,temp)
        next_dist = n^2
        # print("current vertex index:")
        # print(current)
        for v in neighbours:
            alt = current_dist + euclidean_dist(current, v, temp, c, a)
            if alt < dist[v[0],v[1]]:
                dist[v[0],v[1]] = alt
                prev[v[0],v[1]] = current

        # select current with the least value of dist[current]
        flat_min = np.argmin(dist)
        row = flat_min // n
        column = flat_min - row*n
        current = [row, column]

    return prev



def tor2cart(theta, phi, c, a):
    x = (c + a*np.cos(theta)) * np.cos(phi)
    y = (c + a*np.cos(theta)) * np.sin(phi)
    z = a * np.sin(theta)
    return x, y, z

def euclidean_dist(u,v,temp,c,a):
    x1, y1, z1 = tor2cart(temp[u[0]],temp[u[1]],c,a)
    x2, y2, z2 = tor2cart(temp[v[0]],temp[v[1]],c,a)

    return np.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)


def test_find_neighbours():
    n = 100
    c, a = 2, 1
    temp = np.linspace(0, 2*np.pi, n)
    phi, theta = np.meshgrid(temp,temp)
    x,y,z = tor2cart(theta,phi,c,a)

    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, rstride=5, cstride=5, edgecolors='w', alpha = 0.3)
    ax.view_init(36, 26)
    
    # draw vertex
    vertex = [50,0]
    xv,yv,zv = tor2cart(temp[vertex[0]],temp[vertex[1]],c,a)
    ax.scatter(xv, yv, zv, c='r')

    # draw neighbours
    visited = np.zeros((n,n))
    nb_th, nb_ph, nb = find_neighbours(vertex,visited,n,temp)
    x_nb, y_nb,z_nb = tor2cart(nb_th, nb_ph,c,a)
    ax.scatter(x_nb,y_nb,z_nb,c='b')

    plt.show()

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

if __name__ == "__main__":
    n = 100
    c, a = 2, 1
    temp = np.linspace(0, 2*np.pi, n)
    phi, theta = np.meshgrid(temp,temp)
    x,y,z = tor2cart(theta,phi,c,a)

    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, rstride=5, cstride=5, edgecolors='none', alpha = 0.3)
    ax.view_init(36, 26)
    
    # draw two vertices
    start = (79,30)
    end = (0,0)
    xv,yv,zv = tor2cart(temp[np.array([start[0],end[0]])],temp[np.array([start[1],end[1]])],c,a)
    ax.scatter(xv, yv, zv, c='r')

    # get initial path
    prev = Dijkstra(n,start,end,temp)
    current = end
    th_g = [temp[end[0]]]
    ph_g = [temp[end[1]]]
    while not np.array_equal(current,start):
        th_g.append(temp[current[0]])
        ph_g.append(temp[current[1]])
        current = prev[current[0],current[1]]

    # smooth path with fixed point iteration    
    print(th_g)
    print(ph_g)
    th_g = np.asarray(th_g)
    ph_g = np.asarray(ph_g)
    for i in range(10000):
        th_g, ph_g = functional_iteration(th_g,ph_g,a,c)

    xg,yg,zg = tor2cart(th_g,ph_g,c,a)
    ax.plot(xg,yg,zg,c='g')
    plt.show()

    print(th_g)
    print(ph_g)