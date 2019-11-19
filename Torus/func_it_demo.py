import matplotlib.pyplot as plt
from torusmesh import TorusMesh


# functional iteration on path, calculate curve length every iteration

if __name__ == "__main__":

	n = 30
	c, a = 2, 1

	torus = TorusMesh(n,c,a)

	# indices of starting and end point, indices between 0 and n-1
	start = (0,0)
	end = (14,0)

	# draw torus
	fig, ax = torus.draw_torus()

	paths = torus.initial_path(start,end)

	# num_paths = len(paths)
	# prev_paths = []

	# for i in range(num_paths):
	# 	temp, = paths[i].draw_path(ax)
	# 	prev_paths.append(temp)

	# for j in range(200):
	# 	for i in range(num_paths):
	# 		paths[i].functional_iteration()
	# 		prev_paths[i], = paths[i].draw_path(ax, prev_paths[i])
			
	# 		print("curve " + str(i) + " has length " + str(paths[i].curve_length()))
	# 	plt.pause(0.0001)


	p = paths[0]

	prev_path, = p.draw_path(ax)

	curvelength = p.curve_length()
	count = 0
	norm = 10000

	while norm > 10 ** (-5) and count < 1000:
	# while count < 1000:
		p.functional_iteration()
		prev_path, = p.draw_path(ax, prev_path)

		new_curve_length = p.curve_length()
		norm = abs(curvelength - new_curve_length)/curvelength

		curvelength = new_curve_length

		count += 1
		plt.pause(0.0001)
		print("count = " + str(count) + ", length = " + str(curvelength))
