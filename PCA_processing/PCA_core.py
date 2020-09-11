import numpy as np
from scipy.sparse.linalg import eigsh ## no es necesario para la version lite
from scipy.linalg import eigh


# Routines ------------------------------------------------------

def PC_decomposition(data, normal, lite_version = True):  ## sugiero usar un nombre mas apropiado, que especifique lo que hace...

    covariance = np.dot(data.T, data)/np.alen(data)
    #eigenvalues, eigenvectors = eigsh(covariance, k = 100) ## no es necesario para la version lite
    if(lite_version):
        eigenvalues, eigenvectors = eigh(covariance)
    else:
        eigenvalues, eigenvectors = eigsh(covariance, k = 100)
    eigenvectors = eigenvectors[:, np.argsort(-np.abs(eigenvalues))]
    eigenvalues = eigenvalues[np.argsort(-np.abs(eigenvalues))]
    eigenvalues_normalized = eigenvalues/eigenvalues.sum()
    projection = np.dot(eigenvectors.T,data.T)

    return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])

def calc_lite(data, normal):
    ## diferente version del algoritmo
    return 0

## Otras rutinas...

# ------------------------------------------

if __name__ == "__main__":
    print("This file is a collection of routines for PCA, it's not mean to be run isolated.")

