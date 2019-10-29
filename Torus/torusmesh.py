import scipy.io
import numpy as np
from funcpath import FuncPath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go

class TorusMesh:
	"""

	Attributes
	----------
	n: int
		fine-ness of mesh
	half: int
		half of n
	linear: numpy.ndarray
		linear mesh point values
	
	"""

	def __init__(self, n, c, a):
		self.n = n
		self.c = c
		self.a = a

		self.half = int(self.n/2)

		temp = np.linspace(0, 2*np.pi, n+1)
		self.linear = temp[:-1]

		[self.u, self.v] = np.meshgrid(temp,temp)

	
	def tor2cart(self):
		x = (self.c + self.a * np.cos(self.v)) * np.cos(self.u)
		y = (self.c + self.a * np.cos(self.v)) * np.sin(self.u)
		z = self.a * np.sin(self.v)
		return x, y, z

	def shortest_segment(self,a,b): # return shortest path between a and b with mod n
	    if a==b:
	        return []
	    s = np.sign(b-a)
	    temp = np.arange(a,b,s)
	    if temp.size > self.half:
	        return list(np.arange(a,b-s*self.n,-s) % self.n)
	    else:
	        return list(temp % self.n)


	def both_segments(self,a,b): # what does this do again?
	    if a==b:
	        return [[],[]]
	    s = np.sign(b-a)
	    temp = np.arange(a,b,s)

	    return [list(temp), list(np.arange(a,b-s*self.n,-s) % self.n)]

	def initial_path(self, start, end): # generate initial path from two points on the mesh

	    v_a = start[0]
	    v_b = end[0]
	    u_a = start[1]
	    u_b = end[1]
	    

	    # two segment
	    u_seg= self.shortest_segment(u_a,u_b)
	    v_seg1, v_seg2 = self.both_segments(v_a,v_b)

	    if abs(v_a - self.half) < abs(v_b - self.half):
	        # th_a closer to pi
	        v1 = [v_a for i in u_seg] + v_seg1
	        u1 = u_seg + [u_b for i in v_seg1]

	        v2 = [v_a for i in u_seg] + v_seg2
	        u2 = u_seg + [u_b for i in v_seg2]
	    else:
	        # th_b closer to pi
	        v1 = v_seg1 + [v_b for i in u_seg]
	        u1 = [u_a for i in v_seg1] + u_seg

	        v2 = v_seg2 + [v_b for i in u_seg]
	        u2 = [u_a for i in v_seg2] + u_seg

	    # add end point
	    v1 = v1 + [v_b]
	    u1 = u1 + [u_b]
	    v2 = v2 + [v_b]
	    u2 = u2 + [u_b]

	    # three segments: crossing theta = pi
	    v_seg1 = self.shortest_segment(v_a, self.half)
	    v_seg2 = self.shortest_segment(self.half, v_b)

	    v3 = v_seg1 + [self.half for i in u_seg] + v_seg2 + [v_b]
	    u3 = [u_a for i in v_seg1] + u_seg + [u_b for i in v_seg2] + [u_b]

	    path1 = FuncPath(self.linear[np.array(u1)], self.linear[np.array(v1)], self.c, self.a)
	    path2 = FuncPath(self.linear[np.array(u2)], self.linear[np.array(v2)], self.c, self.a)
	    path3 = FuncPath(self.linear[np.array(u3)], self.linear[np.array(v3)], self.c, self.a)
	    
	    return [path1, path2, path3]


	def draw_torus(self):
		# draw the mesh
		x, y, z = self.tor2cart()

		fig = plt.figure()
		ax = fig.gca(projection = '3d')
		ax.set_zlim(-3,3)
		ax.plot_surface(x, y, z, rstride=1, cstride=1, edgecolors='none', alpha = 0.3)
		ax.view_init(36, 26)

		return fig, ax

	def save_data(self):
		x, y, z = self.tor2cart()
		matfile = 'torus_data.mat'
		scipy.io.savemat(matfile, mdict={'xyz': (x,y,z)})

	def plotly_figure(self, filename):
		matdata = scipy.io.loadmat(filename)
		x,y,z =  matdata['xyz']
		return go.Figure(data=[go.Surface(x=x, y=y, z=z, opacity=0.50)])

	def plotly_figure_realtime(self):
		x, y, z = self.tor2cart()
		return go.Figure(data=[go.Surface(x=x, y=y, z=z, opacity=0.50)])