from numpy.linalg import eigh
import numpy as np

def compute_Lsym(A):
    n = A.shape[0]
    m = A.shape[1]
    Anorm = np.zeros((n,n))
    d = np.sum(A,axis = 1)
    dinv = 1./ np.sqrt(d)
    for i in range(n):
        for j in range(m):
            Anorm[i,j] = A[i,j] * dinv[i] * dinv[j]

    I = np.identity(n)
    Lsym = I - Anorm
    return Lsym

def get_symbiharmonic_coords(W):
    d = np.sum(W,axis = 1)
    dinv = 1./ np.sqrt(d)
    Lsym = compute_Lsym(W)
    eigenval, eigenvec = eigh(Lsym)
    inds = np.argsort(eigenval)
    eigenval = eigenval[inds]
    eigenvec = eigenvec[:,inds]
    vol = np.sum(W)
   
   
    inveig = np.square(1./ np.sqrt(eigenval[1:]))
   
   
    Z = np.sqrt(vol) * np.diag(dinv) @ eigenvec[:,1:] @ np.diag(inveig)
    return Z,eigenval
