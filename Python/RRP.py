import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

def RRP(D1,D2):
    n = len(D1)
    r1 = dist_to_rank(D1)
    r2 = dist_to_rank(D2)
    score = np.mean(np.mean(np.absolute(r1 - r2),axis=1))/(n-1)
    x = r1/(n-1)
    y = (r2 - r1)/(n-1)
    result = np.histogram2d(x.flatten(),y.flatten(),bins = [50,100])
    data = gaussian_filter(result[0],sigma=.2)
    X = result[1]
    Y = result[2]
    plt.pcolormesh(X[:-1],Y[:-1],data.T, cmap='inferno', shading='gouraud')
    plt.title('MARE = ' + str(round(score,4)))
    plt.show()
    return score

def dist_to_rank(D):
    n = len(D)
    ind = np.argsort(D,axis=0)
    pointwise_ranks = np.zeros((n-1,n),dtype='int')
    for i in range(n):
        r = np.array(range(1,n+1))
        tmp = np.array(range(1,n+1))
        r[(ind[:,i])] = tmp
        r = r[(setdiff(range(n),[i]))] - 1
        pointwise_ranks[:,i] = r
    return pointwise_ranks
    
def setdiff(lst1,lst2):
    lst3 = [value for value in lst1 if value not in lst2]
    return lst3
