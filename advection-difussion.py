#################################
# This code is created for KTH course Nek5000
# Solve 1-D advection-diffusion problem using SEM
# Zhenyang Yuan  20/03/2021
#################################

# import dependent parkages
import numpy as np
from math import pi
from numpy.linalg import norm
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
    nel = 2  #element number
    L = 1       #domain length
    x = np.linspace(0,L,nel)
    #print(x)
    ord = 3 # choose order of the Legendre polinomial
    (xe,w) = Gauss_Legendre_weight(ord)
    phi = construct_weight_function(xe,ord)
    dphi = differential_weight_function(xe,ord)
    M_e = local_mass(xe,ord,w,phi)
    L_e = local_laplacian(xe,ord,w,dphi)
    D_e = local_derivative(xe,ord,w,phi,dphi)
    A = assembly_matrix(xe,nel)

    #print(phi)

    print("Hello World!")


def construct_weight_function(x,ord):

    phi = np.zeros((ord,len(x)))
    phi_pre = np.zeros((5,len(x)))
    phi_pre[0,:] = np.ones(len(x))
    phi_pre[1,:] = x
    phi_pre[2,:] = 0.5 * (3 * x**2 - 1)
    phi_pre[3,:] = 0.5 * (5 * x**3 - 3 * x)
    phi_pre[4,:] = 0.125 * (35 * x**4 - 30 * x**2 + 3)

    for i in range(ord):
        phi[i,:] = phi_pre[i,:]

    print(phi)
    return phi

def differential_weight_function(xe,ord):
    dphi = np.zeros((ord,len(xe)))
    dphi_pre = np.zeros((5,len(xe)))

    dphi_pre[1,:] = np.ones(len(xe))
    dphi_pre[2,:] = 3 * xe
    dphi_pre[3,:] = 0.5 * (15 * xe**2 - 3)
    dphi_pre[4,:] = 0.125 * (140 * xe**3 - 60 * xe)

    for i in range(ord):
        dphi[i,:] = dphi_pre[i,:]

    #print(dphi)
    return dphi


def Gauss_Legendre_weight(ord):
    xe_pre = [0,[-1/np.sqrt(3),1/np.sqrt(3)],[-np.sqrt(3/5),0,np.sqrt(3/5)],[-np.sqrt(3/7-2/7*np.sqrt(6/5)),np.sqrt(3/7-2/7*np.sqrt(6/5)),-np.sqrt(3/7+2/7*np.sqrt(6/5)),np.sqrt(3/7+2/7*np.sqrt(6/5))]]
    w_pre = [2,[1,1],[5/9,8/9,5/9],[(18+np.sqrt(30))/36,(18+np.sqrt(30))/36,(18-np.sqrt(30))/36,(18-np.sqrt(30))/36]]


    xe = xe_pre[ord-1]
    w  = w_pre[ord-1]
    xe = np.array(xe)  #redefine list to array
    w = np.array(w)
    return (xe,w)


def local_mass(xe,ord,w,phi):

    #print(phi,w)
    M_e = np.zeros((ord,ord))
    for i in range(ord):
        for j in range(ord):
            for k in range(len(w)):
                M_e[i,j] += phi[i,k]*phi[j,k]*w[k]

    #print(M_e)
    return M_e

def local_laplacian(xe,ord,w,dphi):

    L_e = np.zeros((ord,ord))
    for i in range(ord):
        for j in range(ord):
            for k in range(len(w)):
                L_e[i,j] += dphi[i,k]*dphi[j,k]*w[k]

    #print(L_e)
    return L_e



def local_derivative(xe,ord,w,phi,dphi):
    D_e = np.zeros((ord,ord))
    for i in range(ord):
        for j in range(ord):
            for k in range(len(w)):
                D_e[i,j] += phi[i,k]*dphi[j,k]*w[k]

    #print(D_e)
    return D_e


def assembly_matrix(xe,nel):
    #print(np.kron(np.ones((nel,1)),np.eye(xe+2)))  #Kronecker product of two arrays.
    a = np.eye(nel)
    b = np.zeros((len(xe)+2,len(xe)+1))
    b[0:len(xe)+1,0:len(xe)+1] = np.eye(len(xe)+1)
    b[len(xe)+1,len(xe)] = 1
    A = np.zeros(((len(xe)+2)*nel,(len(xe)+2)*nel-nel+1))
    #print(b)


    c = np.kron(a,b)
    c = np.delete(c,(len(xe)+2)*nel-2,0)  #delete the last row to contruct A
    A[0,0] = 1

    A[1:(len(xe)+2)*nel,1:(len(xe)+2)*nel-nel+1] = c
    print(A)


if __name__ == "__main__":
    main()
