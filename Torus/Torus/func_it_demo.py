import functions as fc
import numpy as np
import matplotlib.pyplot as plt

# two particles

# generate initial path

# functional iteration on path, calculate curve length every iteration



if __name__ == "__main__":
	n = 30
	c, a = 2, 1

	temp_mesh = np.linspace(0, 2*np.pi, n+1)
	temp = temp_mesh[:-1]


	[u_torus, v_torus] = np.meshgrid(temp_mesh,temp_mesh)
	x_torus,y_torus,z_torus = fc.tor2cart(u_torus,v_torus,c,a)

	u_particles = [0,2]
	v_particles = [0,2]

	# draw particles on mesh
	x_particles,y_particles,z_particles = fc.tor2cart(u_particles,v_particles,c,a)
	fig, ax = fc.draw_with_vertices(x_torus,y_torus,z_torus,x_particles,y_particles,z_particles)

	plt.show()