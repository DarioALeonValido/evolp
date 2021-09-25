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

lite_version = True
Ngen = 1000
Npc = 20 #the maximum possible number of PCs is the number of genes
tissue_id = "KIRC"
#possible targets:
#tissue_id=["PRAD", "KIRC", "LUSC", "LUAD", "UCEC", "KIRP", "BLCA", "COAD", 
#                   "ESCA", "LIHC", "STAD", "THCA", "BRCA", "HNSC", "READ"]

datapath = "../databases_external/TCGA/"
outputpath = "../databases_generated/TCGA_pca/"  # for Full databases


# Functions ------------------------------------------------------

def read_data(datapath,tissue_id):
    wb = xl.open_workbook(datapath+"sample"+tissue_id+".xls")
    ws = wb.sheet_by_index(0)
    print(ws.nrows,"individual files were identified.")

    sample_type = []
    normal = []
    tumor = []
    data = []
    genes_id = []

    #reading first file to evaluate the cost
    name = ws.cell_value(1,0)
    sample_type.append(ws.cell_value(1,3))
    if(not isfile(datapath+"data"+tissue_id+"/"+name)):
        print(tissue_id, "databases not found, please read", "'"+ datapath + "info.txt'")
        sys.exit()
    file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
    lines = file_object.readlines()
    print("Each file contains",len(lines),"genes.")
    print("Loading large databases...")
    if(ws.cell_value(1,3) == 'Solid Tissue Normal'):
        normal.append(0)
    else:
        tumor.append(0)
    data.append([])
    for line in lines:
        data[-1].append(line.split()[1])
        genes_id.append(line.split()[0].split(sep=".")[0])
    file_object.close()

    #reading the rest of the files
    for i in range(2, ws.nrows):
        name = ws.cell_value(i,0)
        sample_type.append(ws.cell_value(i,3))
        file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
        lines = file_object.readlines()
        if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
            normal.append(i-1)
        else:
            tumor.append(i-1)
        data.append([])
        for line in lines:
            data[-1].append(line.split()[1])
        file_object.close()

    data = np.array(data, dtype="f") + 0.1
    ref = gmean(data[normal])
    data = np.log2(data/ref)

    return([data, normal, tumor, np.array(genes_id)])


# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
if(not isfile(datapath+"sample"+tissue_id+".xls")):
    print(tissue_id, "databases not found, please read", "'"+ datapath + "info.txt'")
    sys.exit()
data, normal, tumor, genes_id = read_data(datapath,tissue_id)
print("Data successfully loaded!")
print("normal =", len(normal))
print("tumor =", len(tumor))
print("total =", len(data))

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
if(lite_version):
    eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PCA_core(data)

    radius_pc1_normal, center_pc1_normal = pca.region(projection[0,normal])
    radius_pc1_tumor, center_pc1_tumor = pca.region(projection[0,tumor])
    radius_pc2_normal, center_pc2_normal = pca.region(projection[1,normal])
    radius_pc2_tumor, center_pc2_tumor = pca.region(projection[1,tumor])
    theta = np.linspace(0, 2*np.pi, 100)

else:
    eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PCA_core(data,True,Npc)

index = np.argpartition(-np.abs(eigenvectors[:, 0]), Npc)[:Npc]
components = eigenvectors[index, 0]
print("Done!")


# Output data and figures --------------------------------------------------

print("Exporting data and plots...")
if(lite_version):
    np.savetxt('pc'+tissue_id+'_dat.dat', projection, fmt='%f')
    np.savetxt('evec'+tissue_id+'_dat.dat', eigenvectors.T, fmt='%f')
    np.savetxt('eval-n'+tissue_id+'_dat.dat', eigenvalues_normalized, fmt='%f')
    np.savetxt('eval'+tissue_id+'_dat.dat', eigenvalues, fmt='%f')
    np.savetxt('evec-t'+tissue_id+'_dat.dat', eigenvectors, fmt='%f')
    np.savetxt('ind20'+tissue_id+'_dat.dat', index, fmt='%i')
    np.savetxt('pc20'+tissue_id+'_dat.dat', components, fmt='%f')
    np.savetxt('ngenes_dat.dat', np.transpose([index,genes_id[index],components]), delimiter="\t", fmt="%s")
    np.savetxt('ind-normal'+tissue_id+'_dat.dat', normal, fmt='%i')
    np.savetxt('ind-tumor'+tissue_id+'_dat.dat', tumor, fmt='%i')

    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    
    ax1.scatter(projection[0,normal], -projection[1,normal], c='b',s=15,label="Normal")
    ax1.scatter(projection[0,tumor], -projection[1,tumor], c='r',s=15,label="Tumor")
    ax1.plot(center_pc1_normal+radius_pc1_normal*np.cos(theta),center_pc2_normal+radius_pc2_normal*np.sin(theta),alpha=0.7)
    ax1.plot(center_pc1_tumor+radius_pc1_tumor*np.cos(theta),center_pc2_tumor+radius_pc2_tumor*np.sin(theta),alpha=0.7)
       
    for i in index:
        ax2.plot([i,i],[0, eigenvectors[i,0]], color='k')

    ax1.set_title(tissue_id), ax1.grid(), ax1.set_xlabel('PC1'), ax1.set_ylabel('-PC2'), ax1.legend()
    ax2.set_title(tissue_id), ax2.grid(), ax2.set_xlabel('Genes'), ax2.set_ylabel('PC1')
    fig1.tight_layout() 
    fig2.tight_layout() 
      
    fig1.savefig(tissue_id+'_PC2_vs_PC1_fig.png')
    fig2.savefig(tissue_id + "_PC1_fig.png")
else:
    np.savetxt(outputpath+'pc'+tissue_id+'.xls', projection.T)

print("Done!")
