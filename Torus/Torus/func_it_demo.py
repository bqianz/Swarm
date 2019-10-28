import functions as fc
import numpy as np
import matplotlib.pyplot as plt
from torusmesh import TorusMesh


# functional iteration on path, calculate curve length every iteration

if __name__ == "__main__":

	n = 30
	c, a = 2, 1

	torus = TorusMesh(n,c,a)

	# indices of starting and end point, indices between 0 and n-1
	start = (0,0)
	end = (5,5)

	# draw torus
	fig, ax = torus.draw_torus()

	p1, p2, p3 = torus.initial_path(start,end)

	prev_plot, = p1.draw_path(ax)

	for i in range(200):
		p1.functional_iteration()
		prev_plot, = p1.draw_path(ax, prev_plot)
		plt.pause(0.0001)



	# TODO: learn animation. is remove() optimal? draw multiple lines at once.