import numpy as np
import matplotlib.pyplot as mtp


# Problem: --- Ut + u*Ux = c Uxx + S
# Not in time for now
# https://web.cecs.pdx.edu/~gerry/class/ME448/lecture/pdf/convectionUpwind.pdf

#Hash properties
Nx = 50
x = np.linspace(0,1,Nx + 1)
dx = 1/Nx

#mid of points (i-dx/2, i , i+dx/2)
xmid = np.array([0.5 * (x[i-1] + x[i]) for i in range(1,len(x))])

# Set u - velocity,c - diffusion, s - source, final time, timestep
u = 1
c = 1
s = 3
tfinal = 1
dt = 0.001

#Set initial condition -- set constant
U = np.array([5 for _ in range(len(xmid))])
lenU = len(U)
t = 0

#Some boundary conditions - now constant
a = 2
b = 4

#while(t < tfinal):

    #Boundary
    #U[0] = a
    #U[-1] = b

    #Diffusion term
    #  east = c_e / dx (U_E - U_P)
    #  west = c_w / dx (U_P - U_W)
    
    #Diff_E = np.array([c/dx * (U[i-1] - U[i]) for i in range(1,lenU)])
    #Diff_W = np.array([c/dx * (U[i] - U[i+1]) for i in range(lenU-1)])

    #Source term
    # Sp dx
    #S = np.array(s * dx for _ in range(lenU))

    #Convection term
    # u_e B_e (U_E - U_P) - u_w B_w (U_W - U_P)
    # B_e = 1/2 when uniform mesh <- B-e =  x_P - x_w / (x_P - x_W)
    #Conv = np.array([u*0.5 * (U[i-1] - U[i]) - u*0.5 * (U[i+1] - U[i]) for i in range(1,lenU-1)])

# Simplyfies to -> -a_e U_e + a_p U_p - a_w U_w = b
# a_e = 1 / dx (D_e - u_e B_e)
# a_w = 1 / dx (D_w - u_w B_w)
# a_p = a_e + a_w
# b = Sp
# D_e = c/dx
# u_e = u
# B_e = 1/2

# matrix A * U = B

A = np.zeros((lenU, lenU))
A[0][0] = 1
A[lenU-1][lenU-1] = 1

for i in range(1,lenU-1):
    a_w = 1/dx * (c/dx + u*0.5)
    a_e = 1/dx * (c/dx - u*0.5)
    A[i][i-1] = -a_w
    A[i][i] = a_e + a_w
    A[i][i+1] = -a_e


B = np.full(lenU, s)
B[0] = a
B[-1] = b

phi = np.linalg.solve(A,B)


mtp.plot(xmid, phi)
mtp.show()