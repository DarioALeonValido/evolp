# -*- coding: utf-8 -*-
"""
The script file applies a Principal Component Analysis (PCA) technique to the
gene expression data that come from Human Protein Atlas (HPA),
http://www.proteinatlas.org/.

For more details on this file see author(s).

@author: JANC
"""
from functools import reduce
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd
from pathlib import Path
#from PCA_core import PCA_core as PCA
from PCA_core_mem import PCAL as PCA  # Principal Component decomposition
                                        # routine with memory boundaries


datapath = Path("../databases_external/HPA/")
fname_normal = Path('normal_tissue.tsv.zip')
fname_tumor = Path('pathology.tsv.zip')

normal = pd.read_table(datapath / fname_normal)
tumor = pd.read_table(datapath / fname_tumor)

# Eliminating some tissues that have few genes correlated in the tissue data
for t in normal.Tissue.drop_duplicates():
    if normal.Gene[normal.Tissue == t].drop_duplicates().shape[0] < 10000:
        normal = normal.drop(normal[normal.Tissue == t].index)
normal.index = np.arange(normal.shape[0])

# Eliminating cases that do not have associates numerical values in the tumor
# data 
tumor = tumor.iloc()[tumor.get([
                                'High', 'Medium', 'Low', 'Not detected'
                                ]).dropna().index]

# Eliminating the genes that are not associated with all the remaining tissues
# in tissues data
G = tumor.Gene.drop_duplicates().to_numpy()
nc = tumor.Cancer.drop_duplicates().shape[0]
for i, g in enumerate(G):
    dummy = tumor[tumor.Gene == g]
    if dummy.shape[0] < nc:
        tumor = tumor.drop(dummy.index)
tumor.index = np.arange(tumor.shape[0])

# Obtaining all the tissues and types of cancer
cancer = tumor.Cancer.drop_duplicates().to_numpy()
tissue = normal.Tissue.drop_duplicates().dropna().to_numpy()

# Intercecting all the sets of genes of every tissue and each type of cancer
# in both sets of data
gene = []
for c in cancer:
    gene.append(tumor.Gene[tumor.Cancer == c].drop_duplicates().to_numpy())
for t in tissue:
    gene.append(normal.Gene[normal.Tissue == t].drop_duplicates().to_numpy())
gene = reduce(np.intersect1d, gene)

# Translating into numerical values the discretized data
matrix_n = np.zeros((tissue.shape[0], gene.shape[0]))
matrix_t = np.zeros((cancer.shape[0], gene.shape[0]))

v = np.array([1, 0, -1, -1])
values = {  # equivalence for data values
    'High': v[0],
    'Medium': v[1],
    'Low': v[2],
    'Not detected': v[3]
}

for i, t in enumerate(tissue):
    dummy = normal[normal.Tissue == t]
    for j, g in enumerate(gene):
        level = dummy.Level[dummy.Gene == g].to_numpy()
        sum_ = n = 0
        for k, l in enumerate(level):
            val = values.get(l, None)
            if val is not None:
                sum_ += val
                n += 1
        matrix_n[i, j] = sum_/n


dict_cancer = {}
for i in cancer:
    dict_cancer[i] = []

for i in range(0, tumor.shape[0], nc):
    if tumor.Gene.iloc()[i] in gene:
        if tumor.iloc()[i:i + nc].get([
                'High', 'Medium', 'Low', 'Not detected'
                ]).isna().to_numpy().any():
            continue
        for j in range(i, i + nc):
            val_l = tumor.iloc()[j].get([
                    'High', 'Medium', 'Low', 'Not detected']).to_list()
            dict_cancer[tumor.Cancer[j]].append(val_l)


for i, c in enumerate(cancer):
    dict_cancer[c] = np.array(dict_cancer[c])
    matrix_t[i] = np.sum(
        (dict_cancer[c].T/dict_cancer[c].sum(axis=1)).T * v, axis=1)

del dict_cancer

# Calculating reference values and normalizing data
ref = matrix_n.mean(axis=0)
matrix = np.concatenate([matrix_n, matrix_t]) - ref
del matrix_n, matrix_t

# Applying the PCA technique to the data
(eigenvalues, eigenvectors,
 eigenvalues_normalized, projection) = PCA(matrix, 50,
                                           matrix_fname='tri_hpa2_dat.dat')

# Getting the Fig. 2
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

n = tissue.shape[0]
ax.scatter(projection[0, :n], projection[1, :n], projection[2, :n],
           label='Normal')
ax.scatter(projection[0, n:], projection[1, n:], projection[2, n:],
           label='Tumor')
ax.set_xlabel('PC1')
ax.set_ylabel('PC2')
ax.set_zlabel('PC3')
ax.legend()
fig.tight_layout()

fig.savefig('Fig2_norm_fig.png')
fig.savefig('Fig2_norm_fig.pdf')
