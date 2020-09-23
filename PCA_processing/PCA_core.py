"""
This module file contains a collection of 
routines for Principal Component Analysis.

For more details on this file see author(s):
JANC, DALV
"""

import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh


# Routines ------------------------------------------------------

def PC_decomp(data, normal, sparse_data = False):

    covariance = np.dot(data.T, data)/np.alen(data)
    if(sparse_data):
        eigenvalues, eigenvectors = eigsh(covariance, k = 100)
    else:
        eigenvalues, eigenvectors = eigh(covariance)
    eigenvectors = eigenvectors[:, np.argsort(-np.abs(eigenvalues))]
    eigenvalues = eigenvalues[np.argsort(-np.abs(eigenvalues))]
    eigenvalues_normalized = eigenvalues/eigenvalues.sum()
    projection = np.dot(eigenvectors.T,data.T)

    return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])


def PC_decomp_mem(data, normal, sparse_data, max_size):

    # Principal Component decomposition routine with memory boundaries

    #return([eigenvalues, eigenvectors, eigenvalues_normalized, projection])
    return 0

def region(data):
    center = data.sum()/np.alen(data)
    radio = np.sqrt(((data - center)**2).sum()/np.alen(data))
    
    return([radio, center])

# ------------------------------------------

if __name__ == "__main__":
    print("This file is a collection of routines for PCA, it's not meant to be run isolated.")
