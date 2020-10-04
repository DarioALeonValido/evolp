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

datapath = '../databases_external/LTEE/Gene_Expression/'
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
        ara = 10**np.array(ara, dtype="f") +0.00001
        anc=[0,1,2,3,8,9,10,11]
        ref = gmean(ara[anc])
        ara = np.log2(ara/ref)

    return gene,ara


# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
gene,ara = read_expression(datapath,ge_file,16)
print('Number of genes: ', len(gene))
print('Number of samples: ', len(ara))
print("Data successfully loaded!")

print("Computing PCA components...")
eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PCA_core(ara)
radius_pc1_anc, center_pc1_anc = pca.region(projection[0,[0,1,2,3,8,9,10,11]])
radius_pc1_evol, center_pc1_evol = pca.region(projection[0,[4,5,6,7,12,13,14,15]])
radius_pc2_anc, center_pc2_anc = pca.region(projection[1,[0,1,2,3,8,9,10,11]])
radius_pc2_evol, center_pc2_evol = pca.region(projection[1,[4,5,6,7,12,13,14,15]])
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
np.savetxt('ind20_dat.dat', index, fmt='%i')
np.savetxt('pc20_dat.dat', components, fmt='%f')
np.savetxt('ngenes_dat.dat', np.transpose([index,gene[index],components]), delimiter="\t", fmt="%s")

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
    
ax1.scatter(-projection[0,[0,1,2,3]], -projection[1,[0,1,2,3]], c='b',s=25,label="Ara+(anc)")
ax1.scatter(-projection[0,[8,9,10,11]],-projection[1,[8,9,10,11]],s=25,label="Ara-(anc)",facecolors='none',edgecolors='b')
ax1.scatter(-projection[0,[4,5,6,7]],-projection[1,[4,5,6,7]], c='r',s=25,label="Ara+(evol)")
ax1.scatter(-projection[0,[12,13,14,15]],-projection[1,[12,13,14,15]],s=25,label="Ara-(evol)",facecolors='none',edgecolors='r')
ax1.plot(-center_pc1_anc+radius_pc1_anc*np.cos(theta),-center_pc2_anc+radius_pc2_anc*np.sin(theta),alpha=0.7)
ax1.plot(-center_pc1_evol+radius_pc1_evol*np.cos(theta),-center_pc2_evol+radius_pc2_evol*np.sin(theta),alpha=0.7)

for i in index:
    ax2.plot([i,i],[0, eigenvectors[i,0]], color='k')

ax1.set_title('Ara'), ax1.grid(), ax1.set_xlabel('PC1'), ax1.set_ylabel('-PC2'), ax1.legend()
ax2.set_title('Ara'), ax2.grid(), ax2.set_xlabel('Genes'), ax2.set_ylabel('PC1')
fig1.tight_layout() 
fig2.tight_layout() 
      
fig1.savefig('_PC2_vs_PC1_fig.png')
fig2.savefig("_PC1_fig.png")

print("Done!")
