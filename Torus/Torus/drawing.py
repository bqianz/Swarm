import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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

def draw_colour(x,y,z,c_array):
    # draw the mesh
    fig = plt.figure()
    ax = fig.gca(projection = '3d')
    ax.set_zlim(-3,3)
    ax.plot_surface(x, y, z, facecolors = c_array)
    ax.view_init(36, 26)
    
    return fig, ax