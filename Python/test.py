import numpy as np
import matplotlib.pyplot as plt
from SASNE import SASNE
from scipy.spatial.distance import pdist,squareform
from RRP import RRP


X = np.loadtxt("../data/imbalanced_test.txt", dtype='f', delimiter=',')


sasne_out = SASNE(X)
embedding = sasne_out[0] # SASNE embedding
Z = sasne_out[1] # original distance coordinates


D1 = squareform(pdist(embedding)) # distances in LD embedding
D2 = squareform(pdist(Z)) # Original graph distances
plt.subplot(1,2,1)
plt.scatter(embedding[:,0],embedding[:,1])
plt.xlabel('SASNE1')
plt.ylabel('SASNE2')
plt.subplot(1,2,2)
res = RRP(D1,D2) # evaluate the embedding with the RRP


plt.show()
