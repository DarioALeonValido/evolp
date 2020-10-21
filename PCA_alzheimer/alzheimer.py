"""
This script file applies the PCA technique to gen expression data
on brain ... that can be found at
https://aging.brain-map.org/download/index
The procedure and main results are described in paper 
arXiv:1706.09813v3 

For more details on this file see author(s):
DALV, AG
"""

import numpy as np
import csv
from scipy.stats.mstats import gmean
import matplotlib.pyplot as plt
import sys
sys.path.append(r'../')
from PCA_cancer import PCA_core as pca
from os.path import *



# General variables ------------------------------------------------------

lite_version = True
Ngen = 1000
Npc = 20 #the maximum possible number of PCs is the number of genes
sample_id = 'FWM' # ['FWM' 'HIP' 'PCx' 'TCx']

datapath = '../databases_external/Aging_Brain/gene_expression/'
donor_file = 'DonorInformation.csv'
samples_file = 'columns-samples.csv'
fpkm_file = 'fpkm_table_normalized.csv'


# Functions ------------------------------------------------------

def read_donor(datapath,donor_file):
    with open(datapath+donor_file) as dfile:
        d_csv = csv.reader(dfile,delimiter=',')

        ND_id = []
        AD_id = []
        for row in d_csv:
            if(row[10] == 'No Dementia'):
                ND_id.append(row[0])
            elif(row[10] == "Alzheimer's Disease Type"):
                AD_id.append(row[0])

    ND_id = np.array(ND_id)
    AD_id = np.array(AD_id)

    return ND_id,AD_id


def read_samples(datapath,samples_file,sample_id,ND_id,AD_id):
    with open(datapath+samples_file) as sfile:
        s_csv = csv.reader(sfile,delimiter=',')

        index_ND = []
        index_AD = []
        i=0
        for row in s_csv:
            if(row[1] in ND_id):
                if(row[8]==sample_id):
                    index_ND.append(i)
            elif(row[1] in AD_id):
                if(row[8]==sample_id):
                    index_AD.append(i)
            i=i+1

    index_ND = np.array(index_ND, dtype="i")
    index_AD = np.array(index_AD, dtype="i")

    return index_ND,index_AD


def read_expression(datapath,fpkm_file,index_ND,index_AD):
    with open(datapath+fpkm_file) as efile:
        e_csv = csv.reader(efile)
        headers = next(e_csv)

        data = []
        genes_id = []
        for row in e_csv:
            data.append([])
            genes_id.append(row[0])
            for c in np.concatenate((index_ND,index_AD)):
                data[-1].append(row[c])

    data = np.array(data, dtype="f")
    data = data.T + 0.1
    ref = gmean(data[:len(index_ND)])
    data = np.log2(data/ref)

    return data,genes_id


# Reading and processing the data ------------------------------------------

print("Reading files from databases...")
ND_id,AD_id = read_donor(datapath,donor_file)
index_ND,index_AD = read_samples(datapath,samples_file,sample_id,ND_id,AD_id)
data,genes_id = read_expression(datapath,fpkm_file,index_ND,index_AD)
len_ND = len(index_ND)
len_AD = len(index_AD)
print(index_ND)
print(index_AD)
print('Number of '+sample_id+' samples:', len_ND+len_AD)
print('Number of ND cases:', len_ND)
print('Number of AD cases:', len_AD)
print('Number of genes:', len(genes_id))
print("Data successfully loaded!")

# reducing the data by Ngen genes for Lite version
if(lite_version and data.shape[1] > Ngen):
    print("Lite version: reducing the number of genes by", Ngen)
    data = data[:, :Ngen]
    genes_id = genes_id[:Ngen]
elif(len(data) < Ngen):
    lite_version = False
    print("Warning: Imposible reduction, number of genes larger than databases.")
    print("The Full version will be performe instead.")

print("Computing PCA components...")
eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PCA_core(data)
radius_pc1_ND,center_pc1_ND = pca.region(projection[0,range(len_ND)])
radius_pc1_AD,center_pc1_AD = pca.region(projection[0,range(len_ND,len_ND+len_AD)])
radius_pc2_ND,center_pc2_ND = pca.region(projection[1,range(len_ND)])
radius_pc2_AD,center_pc2_AD = pca.region(projection[1,range(len_ND,len_ND+len_AD)])
theta = np.linspace(0, 2*np.pi, 100)

index = np.argpartition(-np.abs(eigenvectors[:, 0]), Npc)[:Npc]
components = eigenvectors[index, 0]
print("Done!")


# Output data and figures --------------------------------------------------

print("Exporting data and plots...")
np.savetxt(sample_id+'_pc_dat.dat', projection, fmt='%f')
np.savetxt(sample_id+'_evec_dat.dat', eigenvectors.T, fmt='%f')
np.savetxt(sample_id+'_eval-n_dat.dat', eigenvalues_normalized, fmt='%f')
np.savetxt(sample_id+'_eval_dat.dat', eigenvalues, fmt='%f')
np.savetxt(sample_id+'_evec-t_dat.dat', eigenvectors, fmt='%f')
np.savetxt(sample_id+'_ind20_dat.dat', index, fmt='%i')
np.savetxt(sample_id+'_pc20_dat.dat', components, fmt='%f')
#np.savetxt(sample_id+'ngenes_dat.dat', np.transpose([index,genes_id[index],components]), delimiter="\t", fmt="%s")

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
    
ax1.scatter(projection[0,range(len_ND)], projection[1,range(len_ND)], c='b',s=25,label="ND "+sample_id)
ax1.scatter(projection[0,range(len_ND,len_ND+len_AD)], projection[1,range(len_ND,len_ND+len_AD)], c='r',s=25,label="AD "+sample_id)
ax1.plot(center_pc1_ND+radius_pc1_ND*np.cos(theta),center_pc2_ND+radius_pc2_ND*np.sin(theta),alpha=0.7)
ax1.plot(center_pc1_AD+radius_pc1_AD*np.cos(theta),center_pc2_AD+radius_pc2_AD*np.sin(theta),alpha=0.7)

for i in index:
    ax2.plot([i,i],[0, eigenvectors[i,0]], color='k')

ax1.set_title('Alzheimer'), ax1.grid(), ax1.set_xlabel('PC1'), ax1.set_ylabel('PC2'), ax1.legend()
ax2.set_title('Alzheimer'), ax2.grid(), ax2.set_xlabel('Genes'), ax2.set_ylabel('PC1')
fig1.tight_layout() 
fig2.tight_layout() 
      
fig1.savefig(sample_id+'_PC2_vs_PC1_fig.png')
fig2.savefig(sample_id+"_PC1_fig.png")

print("Done!")
