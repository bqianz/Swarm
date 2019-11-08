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
	h: float
		grid size according to n
	half: int
		half of n
	linear: numpy.ndarray
		linear mesh point values
	c : float
    a : float
        coeffients of torus parametrization. e.g x = (c + a*cos(v)) * cos(u)
	
	"""

	def __init__(self, n, c, a):
		self.n = n
		self.c = c
		self.a = a
		self.h = 2 * np.pi / n

		self.half = int(n/2)

		temp = np.linspace(0, 2*np.pi, n+1)
		self.linear = temp[:-1]

		[self.u, self.v] = np.meshgrid(temp,temp)


	def tor2cart(self):
		"""Calculate cartesian coordinates from toroidal coordinates

        Returns
        -------
        x, y, z : Array of float
        """

		x = (self.c + self.a * np.cos(self.v)) * np.cos(self.u)
		y = (self.c + self.a * np.cos(self.v)) * np.sin(self.u)
		z = self.a * np.sin(self.v)
		return x, y, z

	def shortest_segment(self,a,b):
		"""Return shortest sequence from a to b with mod n.

		Parameters
		----------
		a: int
		b: int
			integer between [0,n)
		
		Returns
		-------
		list
		"""
		if a==b:
			return []
		s = np.sign(b-a)
		temp = np.arange(a,b,s)
		if temp.size > self.half:
			return list(np.arange(a,b-s*self.n,-s) % self.n)
		else:
			return list(temp % self.n)


	def both_segments(self,a,b):
		"""Return both possible sequences from a to b with mod n.
		For example, if a = 1, b = 3, n = 5, then the function returns
		[[1,2,3],[1,0,4,3]]

		Parameters
		----------
		a: int
		b: int
			integer between [0,n)
		
		Returns
		-------
		array of two lists
		"""
		if a==b:
			return [[],[]]
		s = np.sign(b-a)
		temp = np.arange(a,b,s)

		return [list(temp), list(np.arange(a,b-s*self.n,-s) % self.n)]

	def initial_path(self, start, end): # 
		"""Generate three initial paths between two given points on the mesh.

		Parameters
		----------
		start : tuple of lenght 2, data type int
		end : tuple of lenght 2, data type int
			For example, if integer fine-ness n = 20, then elements of tuple can be integers in [0,20)

		Returns
		-------
		path1, path2, path3: array of FuncPath objects (see funcpath.py)
		"""


		u_a = start[0]
		v_a = start[1]

		u_b = end[0]
		v_b = end[1]

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

	# def store_grid(self):
	# 	matfile = 'grid.mat'
	# 	n = self.n

	# 	d = np.zeros(self.[half, n, n], float)
	# 	paths = np.zeros([self.half,n,n], FuncPaths)

	# 	# start = process_time()
	# 	# generate d
	# 	for i in range(half):
	# 		start = (i,0)
	# 		for j in range(n):
	# 			for k in range(n):
	# 				end = (j,k)
	# 				d[i,j,k] = distance(start,end,n,temp)
	# 			print("i,j = " + str(i) + ", " + str(j))
	# 			# print(process_time())

	# 	scipy.io.savemat(matfile, mdict={'distance': d, 'paths': paths})
	# 	return matfile


	def draw_torus(self):
		"""Draw torus as a surface on new figure

		Returns
		-------
		fig, ax: matplotlib variables
		"""
		x, y, z = self.tor2cart()

		fig = plt.figure()
		ax = fig.gca(projection = '3d')
		ax.set_zlim(-3,3)
		ax.plot_surface(x, y, z, rstride=1, cstride=1, edgecolors='none', alpha = 0.3)
		ax.view_init(36, 26)

		return fig, ax

	def save_data(self):
		"""Saves cartesian data of torus surface into a file 'torus_data.mat',
		under variable name 'xyz', in format of a tuple of three numpy arrays.
		"""
		x, y, z = self.tor2cart()
		matfile = 'torus_data.mat'
		scipy.io.savemat(matfile, mdict={'xyz': (x,y,z)})

	def plotly_draw_from_data(self, filename):
		"""Draw torus surface onto new Dash figure from precalculated data.

		Parameters
		----------
		filename: string
			name of file where torus surface data is stored.
		

		Returns
		-------
		Dash graph_objects figure
		"""
		matdata = scipy.io.loadmat(filename)
		x,y,z =  matdata['xyz']
		return go.Figure(data=[go.Surface(x=x, y=y, z=z, opacity=0.50)])

	def plotly_draw_go(self):
		"""
		Draw torus surface onto new Dash figure, graph_objects style.

		Returns
		-------
		Dash graph_objects figure
		"""
		x, y, z = self.tor2cart()
		return go.Figure(
			data=[go.Surface(x=x, y=y, z=z, opacity=0.50)],
			layout=go.Layout(
				uirevision=1
			)
		)

	def plotly_draw(self):
		"""
		Draw torus surface onto new Dash figure, dictionary style.

		Returns
		-------
		Dash figure
		"""
		x, y, z = self.tor2cart()
		return {
			"data": [{"type": "surface",
					  "x": x,
					  "y": y,
					  "z": z,
					  "opacity": 0.5}],
			"layout": {"uirevision": 1} # so that axis don't change between callbacks
			# uirevision is kind of buggy
		}