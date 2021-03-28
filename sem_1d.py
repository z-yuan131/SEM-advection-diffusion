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
from scipy.sparse import csr_matrix
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

    #Initialize varibles
    (nel,L,x,ord,M_g,L_g,D_g,h) = Initialization()
    #calculate GLL points, weights and Legendre polynomial at GLL points
    (xe,w,L_p,dL_p) = Gauss_Legendre_weight(ord)

    x_total = np.zeros((ord-1)*nel+1)           #total points in the domain
    for i in range(nel):
        a = x[i]
        b = x[i+1]
        for j in range(ord-1):
            x_total[(ord-1)*i+j] = xe[j]*(b-a)/2 + (a+b)/2
    x_total[len(x_total)-1] = x[len(x)-1]


    #trial functions and dirivation of trial functions
    phi,dphi = construct_weight_function(xe,ord,L_p,dL_p)
    #build local matrices
    (M_e,L_e,D_e) = local_matrices(ord,w,phi,dphi,h)
    #global assemble matrix and boundary condition matrix
    #note: bc = 0 dirichilet bc = 1 neumann bc = 2 Periodic
    bc = 2
    A,R = assembly_matrix(xe,nel,bc)
    #assemble global matrices
    M_g,L_g,D_g = global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,nel,A,R)


    #solver
    u,uu,x_total = solver(x_total,L,M_g,L_g,D_g,xe,R)
    #plot function
    plot(x_total,u,uu)




def Initialization():
    nel = 5                                     #element number
    ord = 20                                    #number of the Gauss_Lobatto_Legendre points
    L = 1                                       #domain size
    x = np.linspace(0,L,nel+1)                  #grid points

    h = L/nel                                   #element size, even here but can be changed to uneven(more pratical)

    M_g = csr_matrix(( (ord-1)*nel,(ord-1)*nel ))   #using sparse matrix here for the reason that those matrices are diagonal dominated
    L_g = csr_matrix(( (ord-1)*nel,(ord-1)*nel ))
    D_g = csr_matrix(( (ord-1)*nel,(ord-1)*nel ))




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
            P[k+1,:] = ( (2 * k + 1) * x * P[k,:] - k*P[k-1,:] ) / (k + 1);   #legendre poly


        x = xold - ( x * P[N1-1,:] - P[N-1,:] ) / ( N1 * P[N1-1,:] )

    dP[0,:] = 0
    for k in range(N):
        dP[k+1,:] = k*P[k,:] + x*dP[k,:]    #derivatives of legendre poly, not used in this case

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
        # this code can be found in DEVILLE's book, chapter 2, sem Section
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
                M_e[i,j] += phi[k,i]*phi[k,j]*w[k]*h   #since Je = h, dxe/dx = 1/h, some calculation is predone in the code
                L_e[i,j] += dphi[k,i]*dphi[k,j]*w[k]/h
                D_e[i,j] += phi[k,i]*dphi[k,j]*w[k]


    return M_e,L_e,D_e


def assembly_matrix(xe,nel,bc):
    #print(np.kron(np.ones((nel,1)),np.eye(xe+2)))  #Kronecker product of two arrays.

    #element continuity, not very efficient way to contruct A
    a = np.eye(nel)
    b = np.zeros((len(xe),len(xe)-1))
    b[0:len(xe)-1,0:len(xe)-1] = np.eye(len(xe)-1)
    b[len(xe)-1,len(xe)-2] = 1

    A = np.zeros((len(xe)*nel,len(xe)*nel-nel+1))

    c = np.kron(a,b)
    c = np.delete(c,len(xe)*nel-2,0)  #delete the last row to contruct A

    A[0,0] = 1
    A[1:len(xe)*nel,1:len(xe)*nel-nel+1] = c



    #boundary conditions
    N = len(xe)*nel-nel+1
    #if both bc are homogenous dirichilet
    if bc == 0:
        R = np.zeros(( N-2,N ))
        R[:,1:R.shape[1]-1] = np.eye(R.shape[1]-2)
    #if bc on the right is neumann
    if bc == 1:
        R = np.zeros(( N-1,N ))
        R[:,1:R.shape[1]] = np.eye(R.shape[1]-1)
    # if Periodic bc
    if bc == 2:
        R = np.zeros(( N-1,N ))
        R[:,0:R.shape[1]-1] = np.eye(N-1)
        A[A.shape[0]-1,0] = 1
        A[A.shape[0]-1,A.shape[1]-1] = 0


    return A,R


def global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,nel,A,R):
    M_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),M_e)) , A)  #put each element into right place
    L_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),L_e)) , A)
    D_g = np.matmul(np.matmul(np.transpose(A) , np.kron(np.eye(nel),D_e)) , A)


    M_g = np.matmul(np.matmul(R , M_g) , np.transpose(R))
    L_g = np.matmul(np.matmul(R , L_g) , np.transpose(R))
    D_g = np.matmul(np.matmul(R , D_g) , np.transpose(R))

    return M_g,L_g,D_g


def solver_pre(M_g,L_g,D_g,dt,u,xe,c,nv):

    rhs = M_g - c * D_g * dt - nv * L_g * dt               #advection speed = 1, nv = 0.1



    #lhs = reshape_lhs(M_g,xe)                              #useful when doing LU factorization
    #u = solve_banded((len(xe)-1, len(xe)-1), lhs, rhs)     #solve_banded(l_and_u, A, b), system Ax = b
    lhs = np.zeros((M_g.shape[0],M_g.shape[0]))
    for i in range(M_g.shape[0]):
        lhs[i][i] = 1/M_g[i][i]                             #since M_e is diagonal, LU is not needed

    coeffMat = np.matmul(lhs,rhs)                           #pre-calculate coeffcient Matrix

    return coeffMat


def solver(x_total,L,M_g,L_g,D_g,xe,R):
    u_ini = np.transpose(np.sin(x_total*4*3.14/L))
    u = np.matmul(R,u_ini)                                  #apply boundary condition to u
    x_total = np.matmul(R,x_total)                          #also can be done to x to save effort, mathematically, it is right, phically is wrong
    uu = u

    tfinal = 0.1                                            #final time
    c = 1
    nv = 0.1

    dt = stability(M_g,L_g,D_g,c,nv)                   #related to stability analysis, this dt is optimal

    step = math.floor(tfinal/dt)
    

    coeffMat = solver_pre(M_g,L_g,D_g,dt,u,xe,c,nv)


    for t in range(step):
        u = np.matmul(coeffMat,u)
        print('physical time = '+str(dt+t*dt))

    return u,uu,x_total
    
def stability(M_g,L_g,D_g,c,nv):
    from numpy import linalg as LA                          #complex eigenvalue

    AA = nv * L_g + c * D_g
    M_g_n = M_g * 1                                         #not to overwrite
    for i in range(M_g.shape[0]):
        M_g_n[i,i] = 1 / M_g[i,i]
        
    AA = np.matmul(M_g_n,AA)
    w, v = LA.eig(AA)                                       #eigenvalue of stability matrix
    dt = 1
    for i in range(w.shape[0]):
        dt_n = 1 / (w.real[i] ** 2 + w.imag[i] ** 2)
        if dt_n < dt :
            dt = dt_n                                       #find smallest dt, make sure that dt is in that circle
    dt = np.sqrt(dt)
    
    return dt

def plot(x_total,u,uu):

    plt.figure()
    plt.xlabel('x');plt.ylabel('u')
    plt.plot(x_total,u,'-k',lw=2,label='After t')
    plt.plot(x_total,uu,'--r',lw=2,label='initial condition')
    plt.legend(loc='best')
    plt.grid()
    plt.show()



def reshape_lhs(lhs,xe):    #for Lu factorization use, not used here
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

    return lhs_temp




if __name__ == "__main__":
    main()
