import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from construct_graph import construct_graph
from graph_distance import get_symbiharmonic_coords



def SASNE(data):
    W = construct_graph(data)
    res = get_symbiharmonic_coords(W)
    Z = res[0]
    eigenval = res[1]
    init_Y = 1e-4 * Z[:,[1,2]] * np.sqrt(eigenval[1])
    n = len(W)
    perplexity = 0.9 * n
    embedding = TSNE(n_components=2,method='exact',init=init_Y,perplexity=perplexity).fit_transform(Z)

    return embedding


#X_init = np.array([[0,1],[1,0],[1,1],[4,1]])
#X = np.array([[4, 1, 2], [3, 2, 1], [1, 2, 6], [10, 0, 1]])
X = np.loadtxt("imbalanced_test.txt", dtype='f', delimiter=',')
Y = SASNE(X)
print(Y)
plt.scatter(Y[:,0],Y[:,1])
plt.show()


