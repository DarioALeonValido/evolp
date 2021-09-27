# -*- coding: utf-8 -*-
"""
A Principal Component Analysis (PCA) technique with low RAM requirements.

@author: JANC
"""

import numpy as np
from scipy.linalg import eigh
from os.path import isfile

# Covariance Matrix
# =================


def CovMat(data: np.ndarray, buffer_size: int = 8000,
           matrix_fname: str = "upper_triangular_dat.npy",
           fname_cant: str = 'cant_dat.txt') -> None:
    """
    Save in a binary file the upper triangular of the covariance matrix.

    Parameters
    ----------
    data : array_like
        Matrix to which calculary the matrix of covariance.
    buffer_size : int
        Length of de buffer.
    matrix_fname : str
        Name of the file of the covariance matrix.
    fname_cant : str
        Name of the file that contains the number of acces to the disk.
    """
    file = open(matrix_fname, "wb")
    cant = 0
    m, n = data.shape
    buffer = np.empty(buffer_size, dtype=object) if n > buffer_size else\
        np.empty(n, dtype=object)
    for i in range(n):
        buffer[np.mod(i, buffer_size)] = np.dot(data[:, i], data)[i:]/m
        if(np.mod(i+1, buffer_size) == 0):
            np.save(file, np.concatenate(buffer))
            cant += 1
            if(n - i > buffer_size):
                buffer = np.empty(buffer_size, dtype=object)
            else:
                buffer = np.empty(n - i - 1, dtype=object)
    if(buffer.size != 0):
        np.save(file, np.concatenate(buffer))
        cant += 1
    open(fname_cant, 'w').write(str(cant))
    file.close()


# Multiplying matrix for vector
# =============================


def MultMV(vec: np.ndarray, matrix_fname: str) -> np.ndarray:
    """
    Return the matrix product of a vector with de covariance matrix.

    Parameters
    ----------
    vec : array_like
        vector to be multiplied with the covariance matrix
    matrix_fname : str
        Name of the file of the covariance matrix.

    Returns
    -------
    r : array_like
        A vector with the result of multipliplying ``vec`` with the covariance
        matrix
    """
    r = np.zeros_like(vec)
    file = open(matrix_fname, "rb")
    c = np.alen(vec)
    i = 0
    while(i < c):
        buff = np.load(file)
        length = np.alen(buff)
        cumulative = 0
        while(cumulative < length):
            r[i] += np.dot(buff[cumulative:cumulative + c - i], vec[i:])
            r[i+1:] += vec[i] * buff[cumulative + 1:cumulative + c - i]
            cumulative += c - i
            i += 1
        del buff
    file.close()

    return(r)


# Lanzos
# ======


def check(vec: np.ndarray, w: np.ndarray, i: int, lim: float = -1) -> tuple:
    """
    Check if a vector is orthogonal to a set of orthogonal vectors.

    Parameters
    ----------
    vec : array_like
        Set of orthogonal vectors.
    w : array_like
        Vector to be checked for orthogonality.
    i : int
        Number of orthogonal vectors in ``vec``.
    lim : float
        Lim

    Returns
    -------
    out : list
        A list with two elements, the first a boolean that is true if the
        vector is orthogonal with the set of vectors under some approximation
        and the second is a float that indicates the value of the dot product
        between the vector ``w`` and a vector of the vectors of the set
        ``vec``.
    """
    for j in range(1, i + 2):
        buff = np.abs(np.dot(w, vec[j]))
        buff = buff / np.linalg.norm(w)
        if buff > 1e-7:
            return (False, buff)
    return (True, 0)


def ortogonalize(vec: np.ndarray, w: np.ndarray,
                 i: int, buff: float, dim: int) -> np.ndarray:
    """
    Orthogonalize a vector with a set of vectors.

    Parameters
    ----------
    vec : array_like
        Set of orthogonal vectors.
    w : array_like
        Vector to be orthogonalized.
    i : int
        Number of orthogonal vectors in ``vec``.
    buff : float
        Highest allowed value for dot product between vector and vector set.
    dim : int
         Vector size.

    Returns
    -------
    q : array_like
        A vector orthogonal to the vectors of ``vec``.
    """
    q = w.copy()
    for j in range(i + 2):
        q = q - np.dot(q, vec[j]) * vec[j]
    _, buff = check(vec, q, i)

    n = 0
    while(buff > 1e-7):
        q = np.random.randn(dim)
        for j in range(i + 2):
            q = q - np.dot(q, vec[j]) * vec[j]

        _, buff = check(vec, q, i)

        n += 1
        if n == 100:
            print('Something went wrong')
            exit()

    return q


def lanczos(iter: int, dim: int,
            matrix_fname: str = 'upper_triangular_dat.npy') -> tuple:
    """
    Apply Lanczos algorithm to covariance matrix.

    Parameters
    ----------
    iter : int
        Number of iterations of the Lanczos algorithm. Corresponds to the
        dimension of the resulting matrix.
    dim : int
        Covariance matrix dimension.
    matrix_fname : str
        Name of the binary file containing the upper triangular of the
        covariance matrix.

    Returns
    -------
    out : list
        A list of three elements, the first containing the elements of the main
        diagonal, the second the elements of the secondary diagonal and the
        last the vectors used in the Lanczos algorithm.
    """
    vec = np.zeros((iter + 1, dim))
    vec[1] = np.random.rand(dim)
    vec[1] = vec[1]/np.linalg.norm(vec[1])

    alpha = np.zeros(iter)
    beta = np.zeros(iter)

    for i in range(iter - 1):
        w = MultMV(vec[i + 1], matrix_fname=matrix_fname)
        alpha[i] = np.dot(w, vec[i + 1])
        w = w - alpha[i] * vec[i + 1] - beta[i] * vec[i]
        QOrtonormal, valor = check(vec, w, i)
        if not QOrtonormal:
            w = ortogonalize(vec, w, i, valor, dim)

        beta[i + 1] = np.linalg.norm(w)
        vec[i + 2] = w/beta[i + 1]
        i += 1

    w = MultMV(vec[-1], matrix_fname=matrix_fname)
    alpha[-1] = np.dot(w, vec[-1])

    return (alpha, beta[1:], vec)


def PCAL(data: np.ndarray, iterations: int,
         matrix_fname: str = 'upper_triangular_dat.npy',
         buffer_size: int = 8000) -> tuple:
    """
    Subroutine to make the PCA using Lanczos's algorithm.

    Parameters
    ----------
    data : array_like
        Matrix to be analized with PCA.
    iterations : int
        The number of times to be repeated Lanczos's algorithm.
    matrix_fname : str
        Name of the file of the covariance matrix.
    buffer_size : int
        Size of buffer to save in file of covariance matrix

    Returns
    -------
    out : tuple
        tuple with the eigenvalues, eigenvectors, eigenvalues normalized and
        the projections.
    """
    if not isfile(matrix_fname):
        print('Constructing covariance matrix')
        CovMat(data, matrix_fname=matrix_fname, buffer_size=buffer_size)

    a, b, vectors = lanczos(iterations, data.shape[1],
                            matrix_fname=matrix_fname)

    tri = np.diag(a) + np.diag(b, 1) + np.diag(b, -1)

    eigenvalues, eigenvectors = eigh(tri)
    eigenvectors = eigenvectors[:, np.argsort(-np.abs(eigenvalues))]
    eigenvalues = eigenvalues[np.argsort(-np.abs(eigenvalues))]
    eigenvalues_normalized = eigenvalues/eigenvalues.sum()
    e_vec = np.dot(eigenvectors.T, vectors[1:])
    projection = np.dot(e_vec, data.T)

    return (eigenvalues, e_vec, eigenvalues_normalized, projection)
