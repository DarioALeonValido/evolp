# -*- coding: utf-8 -*-
"""
The script file applies a Principal Component Analysis (PCA) technique to the 
gene expression data that come from The Cancer Genome Atlas (TCGA),
https://www.cancer.gov/tcga.

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
#from PCA_core_mem import PCAL as PCA   # Principal Component decomposition 
                                        # routine with memory boundaries


datapath = Path("../databases_generated/tcga-hpa/meanvalues/")

if not isdir(datapath):
    print('You must first unzip the "meanvalues.tar.xy" file, which are located\
 in the "databases_generated/tcga-hpa/" directory, before running this script')
    exit()

fname_tcga = []
for root, dirs, files in walk(datapath):
    for file in files:
        fname_tcga.append(file)

normal = []
tumor = []
for fname in fname_tcga:
    if fname.startswith('normal'):
        normal.append(np.loadtxt(datapath / fname))
    else:
        tumor.append(np.loadtxt(datapath / fname))
normal = np.array(normal)
tumor = np.array(tumor)

ref = gmean(normal)
data = np.concatenate((normal, tumor))/ref
data = np.log2(data)
n = normal.shape[0]
del normal, tumor, ref

eigenvalues, eigenvectors, eigenvalues_normalized, projection = PCAL(
    data, 20, matrix_fname='tri_tcga_5000_dat.npy', buffer_size=4000)

np.savetxt('eigenvalues_dat.dat', eigenvalues)
np.savetxt('eigenvectors_dat.dat', eigenvectors)
np.savetxt('eigenvalues_normalized_dat.dat', eigenvalues_normalized)
np.savetxt('projection_dat.dat', projection)

fig = plt.figure()
ax = Axes3D(fig)

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
