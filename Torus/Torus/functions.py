import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def tor2cart(u, v, c, a):
    x = (c + a*np.cos(v)) * np.cos(u)
    y = (c + a*np.cos(v)) * np.sin(u)
    z = a * np.sin(v)
    return x, y, z

def shortest_segment(a,b,n): # return shortest path between a and b with mod n
    if a==b:
        return []
    half = int(n/2)
    s = np.sign(b-a)
    temp = np.arange(a,b,s)
    if temp.size > half:
        return list(np.arange(a,b-s*n,-s) % n)
    else:
        return list(temp % n)

def initial_path(start,end,n): # generate initial path from two points
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
    th_seg1 = shortest_segment(th_a, half, n)
    th_seg2 = shortest_segment(half, th_b, n)

    th3 = th_seg1 + [half for i in ph_seg] + th_seg2 + [th_b]
    ph3 = [ph_a for i in th_seg1] + ph_seg + [ph_b for i in th_seg2] + [ph_b]

    return temp[np.array(th1)], temp[np.array(ph1)], temp[np.array(th2)], temp[np.array(ph2)], temp[np.array(th3)], temp[np.array(ph3)]

# def iterate(th,ph):
#     curvelength = curve_length(ph,th)
#     count = 0
#     norm = math.inf
#     while norm > 10 ** (-3) and count < 250:
#         th_new, ph_new = functional_iteration(th,ph,a,c)
#         new_curve_length = curve_length(ph_new,th_new)
#         norm = abs(curvelength - new_curve_length)/curvelength
#         curvelength = new_curve_length
#         th = th_new
#         ph = ph_new
#         count += 1
#     return curvelength

# drawing functions
def draw_with_vertices(x,y,z,xv,yv,zv):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, edgecolors='none', alpha = 0.3)
    ax.view_init(36, 26)
    
    # draw 
    ax.scatter(xv, yv, zv, c='r')
    return fig, ax