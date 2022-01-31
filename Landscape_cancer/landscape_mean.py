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
#from PCA_core import PCA_core as PCA
#from PCA_core_mem import PCAL as PCA  # Principal Component decomposition
                                      # routine with memory boundaries

#Dario ----------------------------
import seaborn as sns
import math
import sys
sys.path.append(r'../')
from PCA_cancer import PCA_core as pca
from Manifolds_normal_cancer import PCA_core_mem as pca_mem

tissues_PCA = ['BLCA','BRCA','COAD','ESCA','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','READ','STAD','THCA','UCEC']

def multi_gaussian_dist(mean,dev,coordinates):
    # this is a multivariate gaussian distribution of uncorrelated coordinates
    f = 1.
    s = 1.
    for i in range(len(mean)):
         f=f*math.exp(-((coordinates[i]-mean[i])/dev[i])**2/2.)
         s=s*dev[i]
    f=f/(2.*math.pi*s)

    return f 

def compute_density(pc_data,dev,ranges,npoints):
    axis_pc1 = np.linspace(ranges[0][0],ranges[0][1],npoints[0])
    axis_pc2 = np.linspace(ranges[1][0],ranges[1][1],npoints[1])
    tot_density = np.zeros((len(axis_pc2),len(axis_pc1)))
    for index in range(len(pc_data)):
        for i in range(len(axis_pc1)):
            for j in range(len(axis_pc2)):
                tot_density[j,i]=tot_density[j,i]+multi_gaussian_dist(pc_data[index],dev,[axis_pc1[i],axis_pc2[j]])
    tot_density = tot_density/len(pc_data)

    return axis_pc1, axis_pc2, tot_density

#----------------------------------


datapath = Path("../databases_generated/TCGA_exp/meanvalues/")

if not isdir(datapath):
    print('You must first unzip the "meanvalues.tar.xy" file, which are\
 located in the "databases_generated/tcga-hpa/" directory, before\
 running this script')
    exit()

# Reading the tcga data
tissues = []
normal = []
tumor = []
for root, dirs, files in walk(datapath):
    for file in files:
        if file.startswith('normal'):
            normal.append(np.loadtxt(datapath / file))
        elif file.startswith('tumor'):
            tumor.append(np.loadtxt(datapath / file))
        tissues.append(file)
print(tissues)
normal = np.array(normal)
tumor = np.array(tumor)

# Calculating reference values, normalizing and applying logarithm to the data 
ref = gmean(normal)
data = np.concatenate((normal, tumor))/ref
data = np.log2(data)

n = normal.shape[0]
del normal, tumor, ref

# Applying the PCA technique to the data
eigenvalues, eigenvectors, eigenvalues_normalized, projection = pca_mem.PCAL(
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
fig.savefig('Fig1_fig.pdf')
#plt.show()

#Dario ----------------------------------------------
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.labelsize'] = 17
plt.rcParams['xtick.labelsize'] = 17
plt.rcParams['ytick.labelsize'] = 17
plt.rcParams['legend.fontsize'] = 17
plt.rcParams['lines.markersize'] = 7.

fig2 = plt.figure()
ax2 = fig2.add_subplot()

ranges=[[-150.,200.],[-170.,190.]]
npoints=[36,37]

dev=[25.,40.]
pc1_normal, pc2_normal, density_normal = compute_density(np.transpose([projection[0, :n],-projection[2, :n]]),dev,ranges,npoints)
dev=[25.,40.]
pc1_tumor, pc2_tumor, density_tumor = compute_density(np.transpose([projection[0, n:], -projection[2, n:]]),dev,ranges,npoints)
grid_pc1, grid_pc2 = np.meshgrid(pc1_tumor, pc2_tumor)
tot_density = (density_tumor - density_normal)#/1.1

colormap='RdBu_r'#'seismic'#'coolwarm'
npoints=30
ax2.contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
ax2.contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)

ax2.scatter(projection[0, :n], -projection[2, :n],
           label='Normal',c='b', alpha=0.7)
ax2.scatter(projection[0, n:], -projection[2, n:],
           label='Tumor',c='r', marker="s", alpha=0.7)
#for i in range(n):
#    ax2.annotate(tissues_PCA[i], projection[0][i], -projection[2][i])
#    ax2.annotate(tissues_PCA[i], projection[0][-i], -projection[2][-i])
ax2.legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)

ax2.set_xlabel('PC1')
ax2.set_ylabel('PC3')
plt.tight_layout()
fig2.savefig('Fig5_fig.pdf')

np.savetxt('mean-Normal_dat.dat', np.transpose([tissues, projection[0, :n],projection[1, :n],-projection[2, :n]]), 
delimiter="\t\t", fmt="%s", header='tissue\tPC1\t\tPC2\t\tPC3')
np.savetxt('mean-Tumor_dat.dat', np.transpose([tissues, projection[0, n:],projection[1, n:],-projection[2, n:]]), 
delimiter="\t\t", fmt="%s", header='tissue\tPC1\t\tPC2\t\tPC3')


