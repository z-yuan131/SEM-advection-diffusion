#################################
# This code is created for KTH course Nek5000
# Solve 1-D advection-diffusion problem using SEM
# Zhenyang Yuan  20/03/2021
#################################

# import dependent parkages
import numpy as np
from math import pi
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
    # Initialization
    nel = 20  #element number
    L = 1       #domain length
    x = np.linspace(0,L,nel+1)  #grid point

    ord = 2 # choose order of the Legendre polinomial

    (xe,w) = Gauss_Legendre_weight(ord)
    x_total = np.zeros(len(xe)*nel)
    M_g = np.zeros((len(xe)*nel,len(xe)*nel))
    L_g = np.zeros((len(xe)*nel,len(xe)*nel))
    D_g = np.zeros((len(xe)*nel,len(xe)*nel))

    for i in range(nel):
        a = x[i]
        b = x[i+1]
        phi = construct_weight_function(xe,ord,a,b)
        dphi = differential_weight_function(xe,ord,a,b)

        (M_e,L_e,D_e) = local_matrices(ord,w,phi,dphi,a,b)
        #print(M_e,L_e,D_e)
        #print(D_e)

        for j in range(len(xe)):
            x_total[len(xe)*i+j] = xe[j]*(b-a)/2 + (a+b)/2


        #A = assembly_matrix(xe,nel)

        M_g,L_g,D_g = global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,i,nel)
        #(M_g,L_g,D_g) = global_matrix(M_e,L_e,D_e)

    #print(D_g)

    u = np.transpose(np.sin(x_total*4*3.14/L))
    uu = u
    dt = 0.00001
    step = 20
    tfinal = 0.1

    (lhs,rhs) = solver(M_g,L_g,D_g,dt,u,xe)
    #print(rhs)

    plt.figure()
    plt.plot(x_total,u,'-k',lw=2,label='After t')
    plt.plot(x_total,uu,'--r',lw=2,label='initial condition')
    for t in range(min(np.int_(tfinal/dt),step)):


        b = np.matmul(rhs,u)
        #print(b)
        u = np.matmul(lhs,b)
        #print(u)
        #print(u)

        plt.xlabel('x');plt.ylabel('u')
        plt.plot(x_total,u,'-k',lw=2,label='After t')
        plt.plot(x_total,uu,'--r',lw=2,label='initial condition')
        plt.legend(loc='best')
        plt.grid()
        plt.show()
        #print(u-uu)





    print("Hello World!")

def Gauss_Legendre_weight(ord):
    xe_pre = [0,[-1/np.sqrt(3),1/np.sqrt(3)],[-np.sqrt(3/5),0,np.sqrt(3/5)],[-np.sqrt(3/7-2/7*np.sqrt(6/5)),np.sqrt(3/7-2/7*np.sqrt(6/5)),-np.sqrt(3/7+2/7*np.sqrt(6/5)),np.sqrt(3/7+2/7*np.sqrt(6/5))]]
    w_pre = [2,[1,1],[5/9,8/9,5/9],[(18+np.sqrt(30))/36,(18+np.sqrt(30))/36,(18-np.sqrt(30))/36,(18-np.sqrt(30))/36]]


    xe = xe_pre[ord-1]
    w  = w_pre[ord-1]
    xe = np.array(xe)  #redefine list to array
    w = np.array(w)
    return (xe,w)

def construct_weight_function(xe,ord,a,b):
    #print(xe)
    #xe = (b-a)/2*xe+(a+b)/2*np.ones(len(xe))         #change of interval (mapping)
    #print(xe)
    phi = np.zeros((ord,len(xe)))
    phi_pre = np.zeros((5,len(xe)))
    phi_pre[0,:] = np.ones(len(xe))
    phi_pre[1,:] = xe
    phi_pre[2,:] = 0.5 * (3 * np.power(xe,2) - 1)
    phi_pre[3,:] = 0.5 * (5 * np.power(xe,3) - 3 * xe)
    phi_pre[4,:] = 0.125 * (35 * np.power(xe,4) - 30 * np.power(xe,2) + 3)

    for i in range(ord):
        phi[i,:] = phi_pre[i,:]

    #print(xe,phi)
    return phi

def differential_weight_function(xe,ord,a,b):
    #xe = (b-a)/2*xe+(a+b)/2*np.ones(len(xe))         #change of interval (mapping)

    dphi = np.zeros((ord,len(xe)))
    dphi_pre = np.zeros((5,len(xe)))

    dphi_pre[1,:] = np.ones(len(xe)) * 2 / (b-a)
    dphi_pre[2,:] = 3 * xe * 2 / (b-a)
    dphi_pre[3,:] = 0.5 * (15 * np.power(xe,2) - 3) * 2 / (b-a)
    dphi_pre[4,:] = 0.125 * (140 * np.power(xe,3) - 60 * xe) * 2 / (b-a)

    for i in range(ord):
        dphi[i,:] = dphi_pre[i,:]


    #print(dphi)
    return dphi





def local_matrices(ord,w,phi,dphi,a,b):

    #print(phi,dphi)
    M_e = np.zeros((ord,ord))
    L_e = np.zeros((ord,ord))
    D_e = np.zeros((ord,ord))

    for i in range(ord):
        for j in range(ord):
            for k in range(len(w)):
                M_e[i,j] += phi[i,k]*phi[j,k]*w[k]*(b-a)/2
                L_e[i,j] += dphi[i,k]*dphi[j,k]*w[k]*(b-a)/2
                D_e[i,j] += phi[i,k]*dphi[j,k]*w[k]*(b-a)/2


    #print(phi,dphi)
    #print(L_e)
    return (M_e,L_e,D_e)


# def assembly_matrix(xe,nel):
#     #print(np.kron(np.ones((nel,1)),np.eye(xe+2)))  #Kronecker product of two arrays.
#     a = np.eye(nel)
#     b = np.zeros((len(xe)+2,len(xe)+1))
#     b[0:len(xe)+1,0:len(xe)+1] = np.eye(len(xe)+1)
#     b[len(xe)+1,len(xe)] = 1
#     A = np.zeros(((len(xe)+2)*nel,(len(xe)+2)*nel-nel+1))
#     #print(len(xe))
#
#
#     c = np.kron(a,b)
#     c = np.delete(c,(len(xe)+2)*nel-2,0)  #delete the last row to contruct A
#     A[0,0] = 1
#
#     A[1:(len(xe)+2)*nel,1:(len(xe)+2)*nel-nel+1] = c
#
#     return A

def global_matrix(M_e,M_g,L_e,L_g,D_e,D_g,i,nel):
    #M_g = np.transpose(A) * M_e * A
    A = np.zeros((nel,nel))
    A[i,i] = 1

    M_g += np.kron(A,M_e)
    L_g += np.kron(A,L_e)
    D_g += np.kron(A,D_e)
    print(D_g)
    return M_g,L_g,D_g


def solver(M_g,L_g,D_g,dt,u,xe):

    rhs = M_g -0*D_g * dt - L_g * dt

    #lhs = reshape_lhs(M_g,xe)
    #u = solve_banded((len(xe)-1, len(xe)-1), lhs, rhs)   #solve_banded(l_and_u, A, b), system Ax = b
    lhs = np.zeros((M_g.shape[0],M_g.shape[0]))
    for i in range(M_g.shape[0]):
        lhs[i][i] = 1/M_g[i][i]

    return lhs,rhs

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
