#################################
# This code is created for KTH course Nek5000
# Solve 1-D advection-diffusion problem using SEM
# Zhenyang Yuan  20/03/2021
#################################

# import dependent parkages
import numpy as np
from math import pi
import math
from numpy.linalg import norm
from scipy.linalg import solve_banded
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 15,
          'legend.loc':'best',
          'figure.figsize': (14,5),
          'lines.markerfacecolor':'none',
         'axes.labelsize': 17,
         'axes.titlesize': 17,
         'xtick.labelsize':15,
         'ytick.labelsize':15,
         'grid.alpha':0.6}
pylab.rcParams.update(params)



# main function
def main():


    (nel,L,x,ord,M_g,L_g,D_g,h) = Initialization()
    (xe,w,L_p,dL_p) = Gauss_Legendre_weight(ord)


    x_total = np.zeros((ord-1)*nel+1)
    for i in range(nel):
        a = x[i]
        b = x[i+1]
        for j in range(ord-1):
            x_total[(ord-1)*i+j] = xe[j]*(b-a)/2 + (a+b)/2
    x_total[len(x_total)-1] = x[len(x)-1]




    phi,dphi = construct_weight_function(xe,ord,L_p,dL_p)
    (M_e,L_e,D_e) = local_matrices(ord,w,phi,dphi,h)
    A,R = assembly_matrix(xe,nel)
    M_g,L_g,D_g = global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,nel,A,R)


    #print(dphi)
    #print(D_e,L_e)
    #
    #
    u,uu,x_total = solver(x_total,L,M_g,L_g,D_g,xe,R)
    plot(x_total,u,uu)




def Initialization():
    nel = 5 #element number
    ord = 25 #number of the Gauss_Lobatto_Legendre points
    L = 1 #domain size
    x = np.linspace(0,L,nel+1)   #grid points

    h = L/nel

    M_g = np.zeros(( (ord-1)*nel,(ord-1)*nel ))
    L_g = np.zeros(( (ord-1)*nel,(ord-1)*nel ))
    D_g = np.zeros(( (ord-1)*nel,(ord-1)*nel ))

    return nel,L,x,ord,M_g,L_g,D_g,h




def Gauss_Legendre_weight(ord):
    #  Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
    #  matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
    #  integration and spectral methods.
    #
    #  Reference on LGL nodes and weights:
    #    C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
    #    in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
    #
    #  Written by Greg von Winckel - 04/17/2004


    #calculate GLL points
    N = ord - 1
    N1 = ord;

    m = np.arange(N1)
    #Use the Chebyshev-Gauss-Lobatto nodes as the first guess
    x = np.transpose(np.cos(pi*m/N))

    # The Legendre Vandermonde Matrix
    P = np.zeros((N1,N1))
    dP = np.zeros((N1,N1))

    # Compute P_(N) using the recursion relation
    # Compute its first and second derivatives and
    # update x using the Newton-Raphson method.

    xold = 2
    eps = 10e-10


    while max(abs(x-xold)) > eps:

        xold = x

        P[0,:] = 1
        P[1,:] = x

        for k in range(1,N):
            P[k+1,:] = ( (2 * k + 1) * x * P[k,:] - k*P[k-1,:] ) / (k + 1);


        x = xold - ( x * P[N1-1,:] - P[N-1,:] ) / ( N1 * P[N1-1,:] )

    dP[0,:] = 0
    for k in range(N):
        dP[k+1,:] = k*P[k,:] + x*dP[k,:]

    w = 2 / (N * N1 * P[N1-1,:]**2 );

    xe = -x
    L_p = P
    dL_p = dP

    return (xe,w,L_p,dL_p)



def construct_weight_function(xe,ord,L_p,dL_p):

    phi = np.zeros(( ord,len(xe) ))
    dphi = np.zeros(( ord,len(xe) ))

    for i in range(phi.shape[0]):
        for j in range(phi.shape[1]):
            phi[i][j] = build_lagrange_poly(xe,i,j,0,L_p)  #Li(xej)
            dphi[i][j] = build_lagrange_poly(xe,i,j,1,L_p)  #L'i(xej)
    return phi,dphi



def build_lagrange_poly(xe,i,j,diff,L_p):
    u = 1
    l = 1
    sum = 0

    if diff == 0:
        for m in range(len(xe)):
            if m != i:
                u = u * (xe[j] - xe[m])
                l = l * (xe[i] - xe[m])
        return u / l
    else:
        # for n in range(len(xe)):
        #     if n != i:
        #         for m in range(len(xe)):
        #             if m != i and m != n:
        #                 u = u * (xe[j] - xe[m])
        #                 l = l * (xe[i] - xe[m])
        #         sum += u / l / (xe[i] - xe[n])
        # return sum
        N = len(xe) - 1

        if i == 0 and j == 0:
            return -(N+1)*N/4
        elif i == N and j == N:
            return (N+1)*N/4
        elif i != j:
            return L_p[N,i]/L_p[N,j]/(xe[i] - xe[j])
        else:
            return 0



def local_matrices(ord,w,phi,dphi,h):

    #print(phi,dphi)
    M_e = np.zeros((ord,ord))
    L_e = np.zeros((ord,ord))
    D_e = np.zeros((ord,ord))

    for i in range(ord):
        for j in range(ord):
            for k in range(len(w)):
                # M_e[i,j] += phi[i,k]*phi[j,k]*w[k]*h
                # L_e[i,j] += dphi[i,k]*dphi[j,k]*w[k]/h
                # D_e[i,j] += phi[i,k]*dphi[j,k]*w[k]
                M_e[i,j] += phi[k,i]*phi[k,j]*w[k]*h
                L_e[i,j] += dphi[k,i]*dphi[k,j]*w[k]/h
                D_e[i,j] += phi[k,i]*dphi[k,j]*w[k]


    return M_e,L_e,D_e


def assembly_matrix(xe,nel):
#     #print(np.kron(np.ones((nel,1)),np.eye(xe+2)))  #Kronecker product of two arrays.

    #element continuity
    a = np.eye(nel)
    b = np.zeros((len(xe),len(xe)-1))
    b[0:len(xe)-1,0:len(xe)-1] = np.eye(len(xe)-1)
    b[len(xe)-1,len(xe)-2] = 1
    A = np.zeros((len(xe)*nel,len(xe)*nel-nel+1))
    #print(len(xe))

    c = np.kron(a,b)
    c = np.delete(c,len(xe)*nel-2,0)  #delete the last row to contruct A
    A[0,0] = 1

    A[1:len(xe)*nel,1:len(xe)*nel-nel+1] = c
    #print(A)


    #homogenous dirichilet boundary conditions
    N = len(xe)*nel-nel+1
    #if both bc are dirichilet
    # R = np.zeros(( N-2,N ))
    # R[:,1:R.shape[1]-1] = np.eye(R.shape[1]-2)
    #if bc on the right is neumann
    # R = np.zeros(( N-1,N ))
    # R[:,1:R.shape[1]] = np.eye(R.shape[1]-1)
    # if Periodic bc
    R = np.zeros(( N-1,N ))
    R[:,0:R.shape[1]-1] = np.eye(N-1)
    A[A.shape[0]-1,0] = 1
    A[A.shape[0]-1,A.shape[1]-1] = 0
    # A = A[:,0:A.shape[1]-1]


    return A,R


def global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,nel,A,R):
    M_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),M_e)) , A)
    L_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),L_e)) , A)
    D_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),D_e)) , A)


    M_g = np.matmul(np.matmul(R , M_g) , np.transpose(R))
    L_g = np.matmul(np.matmul(R , L_g) , np.transpose(R))
    D_g = np.matmul(np.matmul(R , D_g) , np.transpose(R))

    return M_g,L_g,D_g


def solver_pre(M_g,L_g,D_g,dt,u,xe):

    rhs = M_g - 1 * D_g * dt - 0.1 * L_g * dt



    #lhs = reshape_lhs(M_g,xe)   #useful when doing LU factorization
    #u = solve_banded((len(xe)-1, len(xe)-1), lhs, rhs)   #solve_banded(l_and_u, A, b), system Ax = b
    lhs = np.zeros((M_g.shape[0],M_g.shape[0]))
    for i in range(M_g.shape[0]):
        lhs[i][i] = 1/M_g[i][i]

    coeffMat = np.matmul(lhs,rhs)




    return coeffMat


def solver(x_total,L,M_g,L_g,D_g,xe,R):
    u_ini = np.transpose(np.sin(x_total*4*3.14/L))
    u = np.matmul(R,u_ini)
    x_total = np.matmul(R,x_total)
    uu = u

    tfinal = 0.1
    CFL = 0.0001
    c = 1
    dt = CFL*(x_total[2]-x_total[1])/c
    #dt = 0.000001
    step = math.floor(tfinal/dt)
    # step = 1

    coeffMat = solver_pre(M_g,L_g,D_g,dt,u,xe)


    for t in range(step):
        u = np.matmul(coeffMat,u)

    return u,uu,x_total

def plot(x_total,u,uu):

    plt.figure()
    plt.xlabel('x');plt.ylabel('u')
    plt.plot(x_total,u,'-k',lw=2,label='After t')
    plt.plot(x_total,uu,'--r',lw=2,label='initial condition')
    plt.legend(loc='best')
    plt.grid()
    plt.show()



def reshape_lhs(lhs,xe):
    lhs_temp = np.zeros((len(xe)*2-1,lhs.shape[0]))

    for i in range(lhs.shape[0]):
        if i < lhs.shape[0] - len(xe):
            for j in range(len(xe)):
                lhs_temp[len(xe)-1-j,i+j] = lhs[i,i+j]
                lhs_temp[len(xe)-1+j,i] = lhs[i+j,i]
        else:
            for j in range(lhs.shape[0] - i):
                lhs_temp[len(xe)-1-j,i+j] = lhs[i,i+j]
                lhs_temp[len(xe)-1+j,i] = lhs[i+j,i]

#    print(rhs_temp)
    return lhs_temp




if __name__ == "__main__":
    main()
