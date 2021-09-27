# -*- coding: utf-8 -*-
"""
The script file applies a Principal Component Analysis (PCA) technique to the
gene expression data that come from The Cancer Genome Atlas (TCGA),
https://www.cancer.gov/tcga/.

For more details on this file see author(s).

@author: JANC
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from os import walk
from os.path import isdir
from pathlib import Path
from scipy.stats import gmean
from PCA_core import PCA_core as PCA
# from PCA_core_mem import PCAL as PCA  # Principal Component decomposition
                                        # routine with memory boundaries


datapath = Path("../databases_generated/tcga-hpa/meanvalues/")

if not isdir(datapath):
    print('You must first unzip the "meanvalues.tar.xy" file, which are\
 located in the "databases_generated/tcga-hpa/" directory, before\
 running this script')
    exit()

# Reading the tcga data
normal = []
tumor = []
for root, dirs, files in walk(datapath):
    for file in files:
        if file.startswith('normal'):
            normal.append(np.loadtxt(datapath / file))
        elif file.startswith('tumor'):
            tumor.append(np.loadtxt(datapath / file))
normal = np.array(normal)
tumor = np.array(tumor)

# Calculating reference values, normalizing and applying logarithm to the data 
ref = gmean(normal)
data = np.concatenate((normal, tumor))/ref
data = np.log2(data)

n = normal.shape[0]
del normal, tumor, ref

# Applying the PCA technique to the data
eigenvalues, eigenvectors, eigenvalues_normalized, projection = PCA(
    data, 20, matrix_fname='tri_tcga_5000_dat.npy', buffer_size=4000)

# Getting the top panel of Fig. 1
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

ax.scatter(-projection[0, :n], -projection[1, :n], projection[2, :n],
           label='Normal')
ax.scatter(-projection[0, n:], -projection[1, n:], projection[2, n:],
           label='Tumor')
ax.legend()
ax.set_xlabel('- PC1')
ax.set_ylabel('- PC2')
ax.set_zlabel('PC3')
ax.view_init(18, -76)
fig.savefig('Fig1_fig.png')
fig.savefig('Fig1_fig.pdf')
plt.show()
