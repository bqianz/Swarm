

# global variables
a = 1
c = 2
n = 50



# newtons
def scaled_norm(dY): # returned scaled norm of dY
    l = length(dY)
    sum1 = abs(dY[0]) + abs(dY[1])
    sum2 = abs(dY[2]) + abs(dY[3])

    for i in range(4,l,4):
        sum1 += abs(dY[i+2]) + abs(dY[i+3])
        sum2 += abs(dY[i]) + abs(dY[i+1])

    return sum1 + sum2/10

def make_Y_intial(u,v): # returns Y of lenght 4m

    # form Y for from initial paths on the mesh
    m = lenght(u)
    Y = np.zeros(4*m)

    Y[0] = u[0]
    Y[1] = v[0]

    for i in range(1,m):

        # p and q with arclength parametrization
        u_diff = u[i] - u[i-1]
        if(u_diff == 0):
            Y[4*i] = 0
            Y[4*i+1] = np.sign(v[i] - v[i-1]) / a

        else:
            Y[4*i]= = np.sign(u_diff) / c + a * cos(v[i])

        # u and v
        Y[4*i+2] = u[i]
        Y[4*i+3] = v[i]
    return Y

def make_arrays_from_Y(Y): # returns v,p,q
    m = int(length(Y)/4)

    u = np.zeros(m)
    v = np.zeros(m)
    p = np.zeros(m)
    q = np.zeros(m)

    u[0] = Y[0]
    v[0] = Y[1]
    p[0] = Y[2]
    q[0] = Y[3]
    for i in range(1,m):
        u[i] = Y[4*i+2]
        v[i] = Y[4*i+3]
        p[i] = Y[4*i]
        q[i] = Y[4*i+1]

    return u,v,p,q

def make_values(Y): # returns om1,om2,ps1,h
    # generate omega1, omega2, psi1

    m = int(length(Y) / 4)
    u,v,p,q = make_arrays_from_Y(Y)

    # generate array of christoffel symbols, size m
    temp1 = c + a*np.cos(v)

    chr11 = np.sin(v) * temp1 / a
    chr12 = a * np.sin(v) / temp1

    om1 = chr12 * q
    om2 = chr11 * p
    ps1 = chr12 * p

    u_diff = u - np.concatenate(([u[0]],u[1:]))
    v_diff = v - np.concatenate(([v[0]],v[1:]))

    h = np.sqrt( np.square(temp1 * u_diff) + np.square( a * v_diff) )

    return om1, om2, ps1, h

def make_Jacobian(Y): # returns jacobian matrix J
    m4 = length(Y)
    m = int(m4/4)
    # generate array 
    J = np.zeros((m4,m4))

    om1,om2,ps1,h = make_values(Y)

    h2 = - h/2

    # make B_k for k = 2,..,m
    for k in range(1,m):
        pr = k*4-2 # pivot row
        pc = k*4 # pivot column

        h_k = h[k]

        J[pr,pc] = h2[k]
        J[pr+1,pc+1] = h2[k]
        J[pr,pc+2] = 1
        J[pr+1,pc+3] = 1
        J[pr+2,pc] = h_k*om1[k] + 1
        J[pr+2,pc+1] = h_k * ps1[k]
        J[pr+3,pc] = h_k*om2[k]
        J[pr+3,pc+1] = 1

    # make A_k for k=3,...,m
    for k in range(2,m):
        pr = (k-1)*4 + 2
        pc = (k-1)*4

        h_k = h[k]

        J[pr,pc] = h2[k]
        J[pr+1,pc+1] = h2[k]
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

        J[pr,pc+2] = h2[k]
        J[pr+1,pc+3] = h2[k]

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

def determine_mu(s_norm): # returns correction scalar mu

    temp_val = s_norm / l

    if temp_val <= 0.005:
        mu = 0.6
    else if temp_val <= 0.05:
        mu = 0.4
    else:
        mu = 0.2

    return mu

def newtons(Y): # returns new Y
    m4 = length(Y)
    J = make_Jacobian(Y)
    dY = LA.lstsq(J,np.zeros(m4))
    mu = determine_mu(dY)
    return Y + mu * dY

def it_newtons(u,v):
    count = 0
    Y = make_Y_intial(u,v)
    l = 4 * length(u)
    F = np.zeros(l)

    # continue here: specify stopping condition 
    while count < 50:
        J = make_Jacobian(Y)
        dY = LA.lstsq(J,F)
        sn = scaled_norm(dY)
        mu = determine_mu(sn)
        Y += mu * dY
        count += 1
    return curvelength

