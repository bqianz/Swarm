class NewtonPath:
	# attributes: length, m, Y

	def __init__(self,u,v):
		# u, v of length m, Y of length 4m
		self.m = len(u)
		self.length = 4*self.m

		# calculate data
		self.Y = np.zeros(self.length)
	    self.Y[0] = u[0]
	    self.Y[1] = v[0]

	    for i in range(1,m):

	        # p and q with arclength parametrization
	        # why does the parametrization have to be arclength?
	        u_diff = u[i] - u[i-1]
	        if(u_diff == 0):
	            self.Y[4*i] = 0
	            self.Y[4*i+1] = np.sign(v[i] - v[i-1]) / a

	        else:
	            self.Y[4*i]= = np.sign(u_diff) / c + a * cos(v[i])

	        # u and v
	        self.Y[4*i+2] = u[i]
	        self.Y[4*i+3] = v[i]

	def precalculate_for_Jacobian(self): # returns om1,om2,ps1,h

		u = np.concatenate(([self.Y[0]], self.Y[6::4]))
		v = np.concatenate(([self.Y[1]], self.Y[7::4]))
		p = np.concatenate(([self.Y[2]], self.Y[4::4]))
		q = np.concatenate(([self.Y[3]], self.Y[5::4]))

	    # generate array of christoffel symbols, size m
	    temp = c + a*np.cos(v)

	    chr11 = np.sin(v) * temp / a
	    chr12 = a * np.sin(v) / temp

	    om1 = chr12 * q
	    om2 = chr11 * p
	    ps1 = chr12 * p

	    u_diff = u - np.concatenate(([u[0]],u[1:]))
	    v_diff = v - np.concatenate(([v[0]],v[1:]))

	    h = np.sqrt( np.square(temp * u_diff) + np.square( a * v_diff) ) # is this right?

	    return om1, om2, ps1, h

	def make_Jacobian(self): # returns jacobian matrix J

	    J = np.zeros((self.length,self.length))

	    om1,om2,ps1,h = precalculate_for_Jacobian(self.Y)

	    half = - h/2

	    # make B_k for k = 2,..,m
	    for k in range(1,self.m):
	        pr = k*4-2 # pivot row
	        pc = k*4 # pivot column

	        h_k = h[k]

	        J[pr,pc] = half[k]
	        J[pr+1,pc+1] = half[k]
	        J[pr,pc+2] = 1
	        J[pr+1,pc+3] = 1
	        J[pr+2,pc] = h_k*om1[k] + 1
	        J[pr+2,pc+1] = h_k * ps1[k]
	        J[pr+3,pc] = h_k*om2[k]
	        J[pr+3,pc+1] = 1

	    # make A_k for k=3,...,m
	    for k in range(2,self.m):
	        pr = (k-1)*4 + 2
	        pc = (k-1)*4

	        h_k = h[k]

	        J[pr,pc] = half[k]
	        J[pr+1,pc+1] = half[k]
	        J[pr,pc+2] = -1
	        J[pr+1,pc+3] = -1
	        J[pr+2,pc] = h_k*om1[k-1] - 1
	        J[pr+2,pc+1] = h_k * ps1[k-1]
	        J[pr+3,pc] = h_l * om2[k-1]
	        J[pr+3,pc+1] = -1

	    # make A_2
	    for k in range(1,2):
	        pr = 2
	        pc = 0

	        h_k = h[k]

	        J[pr,pc] = -1
	        J[pr+1,pc+1] = -1

	        J[pr,pc+2] = half[k]
	        J[pr+1,pc+3] = half[k]

	        J[pr+2,pc+2] = h_k*om1[k-1] - 1
	        J[pr+2,pc+3] = h_k * ps1[k-1]
	        J[pr+3,pc+2] = h_l * om2[k-1]
	        J[pr+3,pc+3] = -1

	    # make B_1
	    J[0,0] = 1
	    J[1,1] = 1

	    # make A_m+1
	    J[m-1,m-1] = 1
	    J[m-2,m-2] = 1

	    return J

	def determine_mu(self,dY): # returns correction scalar mu

	    sum1 = abs(dY[0]) + abs(dY[1])
	    sum2 = abs(dY[2]) + abs(dY[3])

	    for i in range(1,self.m):
	        sum1 += abs(dY[4*i+2]) + abs(dY[4*i+3])
	        sum2 += abs(dY[4*i]) + abs(dY[4*i+1])

	    s_norm = sum1 + sum2/10

	    temp_val = s_norm / self.length

	    if temp_val <= 0.005:
	        mu = 0.6
	    else if temp_val <= 0.05:
	        mu = 0.4
	    else:
	        mu = 0.2

	    return mu

	def newtons(self): # updates new self.Y
	    J = self.make_Jacobian()
	    dY = LA.lstsq(J,np.zeros(self.length))
	    mu = self.determine_mu(dY)
	    self.Y = self.Y + mu * dY

	# def does_it_converge(self):

