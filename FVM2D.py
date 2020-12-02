import numpy as np
import matplotlib.pyplot as mpt

# Problem: --- Ut + u*Ux = c Uxx + S
# Not in time for now

def createMesh2D(Nx, Ny, xlen, ylen):
    """
    Input:
    -Nx = number of cells in x
    -Ny = number of cells in y
    -xlen = lenght in x (whole object)
    -ylen = lenght in y (whole object)

    Output:
    -X,Y = node locations, even boundaries
    -XW,YW = interface locations on east/south of each node
    """

    dx = xlen/Nx
    dy = ylen/Ny
    XW = np.zeros(Nx+2)
    for i in range(2,Nx+2):
        XW[i] = XW[i-1] + dx
    YW = np.zeros(Ny+2)
    for i in range(2,Ny+2):
        YW[i] = YW[i-1] + dy

    X = np.zeros(len(XW))
    for i in range(1,Nx+1):
        X[i] = (XW[i] + XW[i+1]) / 2
    X[Nx+1] = xlen

    Y = np.zeros(len(YW))
    for i in range(1,Ny+1):
        Y[i] = (YW[i] + YW[i+1]) / 2
    Y[Ny+1] = ylen

    print(XW)
    print(YW)
    print(X)
    print(Y)

def centralDiff2D(parameter_list):
    """
    Input:
    -u = uniform velocity (scalar)
    -gam = uniform diffusion coefficient (scalar)
    -X,Y = cell centers
    -XW,YW = positions of faces west/north
    -phib = boundary values
    -src = source term
    
    Output:
    -A = matrix
    -b = right hand side
    """
    pass



createMesh2D(5,5,1,1)

"""

Nx = 5
X,Y = np.mgrid[-5:5:complex(0,Nx), -5:5:5j]
Nodes = np.zeros((Nx, Nx))
tmp = 0
for i in range(Nx):
    for j in range(Nx):
        Nodes[i][j] = tmp
        tmp += 1

dx = 1/Nx
lenU = Nx * Nx
# Set u - velocity,gamma - diffusion, s - source, final time, timestep
vel = 1
gamma = 1
source = 3

#Some boundary conditions - now constant
a = 2 #RIGHT
b = 4 #TOP
c = 3 #LEFT
d = 5 #BOTTOM

#------JUST DIFFUSION---Might be little bit wrong -- needs change, https://en.wikipedia.org/wiki/Finite_volume_method_for_two_dimensional_diffusion_problem

A = np.zeros((lenU, lenU))
borders = []

#TODO: this shit - needs border nodes from nodes
for i in range(Nx):
    for j in range(Nx):
        if(i == 0 or j == 0 or i == Nx-1 or j == Nx-1):
            borders.append(Nodes[i][j])

#print(Nodes)
#print(borders)

for i in range(0,lenU):
    if(i in borders):
        A[i][i] = 1
        continue

    a_w = gamma * dx / (dx/2)
    a_e = gamma * dx / (dx/2)
    a_s = gamma * dx / (dx/2)
    a_n = gamma * dx / (dx/2)
    a_p = a_w + a_e + a_s + a_n - source

    A[i][i-1] = -a_w
    A[i][i] = a_p
    A[i][i+1] = -a_e
    A[i+1][i] = -a_s
    A[i-1][i] = -a_n


B = np.full(lenU, 1)
for i in borders:
    if(i < Nx):
        B[int(i)] = d
    elif(i >= lenU - Nx):
        B[int(i)] = b
    else:
        B[int(i)] = c

print(A)
print(B)

U = np.linalg.solve(A,B)

print(U)
print(U.shape)
U = np.reshape(U,(Nx,Nx))

print(U)

"""