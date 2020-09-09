import numpy as np
from scipy.stats.mstats import gmean
#from scipy.sparse.linalg import eigsh ## no es necesario para la version lite
from scipy.linalg import eigh


# Routines ------------------------------------------------------

def calc(data, normal):  ## sugiero usar un nombre mas apropiado, que especifique lo que hace...
    #ref = gmean(data[normal])  ## esto no me funciono asi
    data_n = []
    for i in normal:
        data_n.append(data[i])
    ref = gmean(data_n)

    t = np.log2(data/ref)
    covariance = np.dot(t.T, t)/t.shape[0]
    #eigenvalues, eigenvectors = eigsh(covariance, k = 100) ## no es necesario para la version lite
    eigenvalues, eigenvectors = eigh(covariance)
    eigenvalues_normalized = eigenvalues/eigenvalues.sum()
    projection = np.dot(eigenvalues.T,t.T)
    
    return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])

def calc_lite(data, normal):
    ## diferente version del algoritmo  
    return 0

## Otras rutinas...

# ------------------------------------------

if __name__ == "__main__":
    print("This file is a collection of routines for PCA, it's not mean to be run isolated.")

