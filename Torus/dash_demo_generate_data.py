from torusmesh import TorusMesh
import numpy as np
import scipy.io

n = 50
c, a = 2, 1

torus = TorusMesh(n,c,a)
torus_data = torus.tor2cart()

start = [0,0]
end = [8,8]
path, _, _ = torus.initial_path(start,end)

num_iterations = 60

x,y,z = path.tor2cart()
x_edge = x[[0,-1]]
y_edge = y[[0,-1]]
z_edge = z[[0,-1]]
# data[0] contains u, data[1]  contains v 

data = np.zeros([3, num_iterations + 1, path.length])
data[0,0], data[1,0], data[2,0] = path.tor2cart()

for i in range(1, num_iterations+1):
    path.functional_iteration()
    data[0,i], data[1,i], data[2,i] = path.tor2cart()

matfile = 'dash_demo_data.mat'
scipy.io.savemat(matfile, mdict={
    'path_data': data,
    'torus_data': torus_data,
    'n' : n,
    'c' : c,
    'a' : a,
    'end_points': [x_edge, y_edge, z_edge]
    })


