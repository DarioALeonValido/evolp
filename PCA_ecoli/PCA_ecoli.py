"""
This script file applies the PCA technique to the TCGA data at 
https://www.cancer.gov/tcga
The procedure and main results are described in paper 
arXiv:1706.09813v3 
A more general analysis can be found in paper
arXiv:2003.07828v3 

For more details on this file see author(s):
JANC, DALV
"""

import numpy as np
import xlrd as xl
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
import sys
sys.path.append(r'../')
from PCA_cancer import PCA_core as pca
from os.path import *



# General variables ------------------------------------------------------

Npc = 20 #the maximum possible number of PCs is the number of genes

datapath = '../databases_external/LTEE/genes/'
ge_file = 'expression20k.dat'


# Functions ------------------------------------------------------

def read_expression(datapath,ge_file,ncol):
    with open(datapath+ge_file, 'r') as gefile:
        lines = gefile.readlines()

        gene = []
        ara = []
        for j in range(ncol):
            ara.append([])   
        for line in lines[1:]:
            gene.append(line.split()[0])
            for j in range(ncol):
                ara[j].append(line.split()[j+1])     
        gene = np.array(gene)
        ara = -np.array(ara, dtype="f") +0.0001
        araP = ara[:8]
        anc=[0,1,2,3,8,9,10,11]
        ref = gmean(ara[anc])
        ara = np.log2(ara/ref)
        ancP=[0,1,2,3]
        refP = gmean(araP[ancP])
        araP = np.log2(araP/refP)

    return gene,ara,araP


# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
gene,ara,araP = read_expression(datapath,ge_file,16)
print('Number of genes: ', len(gene))
print('Number of samples: ', len(ara))
print("Data successfully loaded!")

print("Computing PCA components...")
eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PC_decomp(ara)
eigenvaluesP,eigenvectorsP,eigenvalues_normalizedP,projectionP = pca.PC_decomp(araP)
radius_pc1_normal, center_pc1_normal = pca.region(projection[0,[0,1,2,3,8,9,10,11]])
radius_pc1_tumor, center_pc1_tumor = pca.region(projection[0,[4,5,6,7,12,13,14,15]])
radius_pc2_normal, center_pc2_normal = pca.region(projection[1,[0,1,2,3,8,9,10,11]])
radius_pc2_tumor, center_pc2_tumor = pca.region(projection[1,[4,5,6,7,12,13,14,15]])
theta = np.linspace(0, 2*np.pi, 100)

index = np.argpartition(-np.abs(eigenvectors[:, 0]), Npc)[:Npc]
components = eigenvectors[index, 0]
print("Done!")


# Output data and figures --------------------------------------------------

print("Exporting data and plots...")
np.savetxt('pc_dat.dat', projection, fmt='%f')
np.savetxt('evec_dat.dat', eigenvectors.T, fmt='%f')
np.savetxt('eval-n_dat.dat', eigenvalues_normalized, fmt='%f')
np.savetxt('eval_dat.dat', eigenvalues, fmt='%f')
np.savetxt('evec-t_dat.dat', eigenvectors, fmt='%f')
#np.savetxt('ind20'+tissue_id+'_dat.dat', index, fmt='%i')
#np.savetxt('pc20'+tissue_id+'_dat.dat', components, fmt='%f')
np.savetxt('ngenes_dat.dat', np.transpose([index,gene[index],components]), delimiter="\t", fmt="%s")
#np.savetxt('ind-normal'+tissue_id+'_dat.dat', normal, fmt='%i')
#np.savetxt('ind-tumor'+tissue_id+'_dat.dat', tumor, fmt='%i')

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig1P, ax1P = plt.subplots()
    
ax1.scatter(projection[0,[0,1,2,3]], -projection[1,[0,1,2,3]], c='b',s=25,label="Ara+(anc)")
ax1.scatter(projection[0,[8,9,10,11]],-projection[1,[8,9,10,11]],s=25,label="Ara-(anc)",facecolors='none',edgecolors='b')
ax1.scatter(projection[0,[4,5,6,7]],-projection[1,[4,5,6,7]], c='r',s=25,label="Ara+(evol)")
ax1.scatter(projection[0,[12,13,14,15]],-projection[1,[12,13,14,15]],s=25,label="Ara-(evol)",facecolors='none',edgecolors='r')
ax1.plot(center_pc1_normal+radius_pc1_normal*np.cos(theta),center_pc2_normal+radius_pc2_normal*np.sin(theta),alpha=0.7)
ax1.plot(center_pc1_tumor+radius_pc1_tumor*np.cos(theta),center_pc2_tumor+radius_pc2_tumor*np.sin(theta),alpha=0.7)

ax1P.scatter(projectionP[0,[0,1,2,3]], projectionP[1,[0,1,2,3]], c='b',s=15,label="Ara+(anc)")
ax1P.scatter(projectionP[0,[4,5,6,7]], -projectionP[1,[4,5,6,7]], c='r',s=15,label="Ara+(evol)")

for i in index:
    ax2.plot([i,i],[0, eigenvectors[i,0]], color='k')

ax1.set_title('Ara'), ax1.grid(), ax1.set_xlabel('PC1'), ax1.set_ylabel('-PC2'), ax1.legend()
ax2.set_title('Ara'), ax2.grid(), ax2.set_xlabel('Genes'), ax2.set_ylabel('PC1')
fig1.tight_layout() 
fig2.tight_layout() 
      
fig1.savefig('_PC2_vs_PC1_fig.png')
fig2.savefig("_PC1_fig.png")
fig1P.savefig('AraP_PC2_vs_PC1_fig.png')

print("Done!")
