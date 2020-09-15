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
from pathlib import Path
import shutil
import PCA_core as pca


# Names ----------------------------------------------------------

lite_version = True
Ngen = 1000
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

    #reading first file to evaluate the cost
    name = ws.cell_value(1,0)
    sample_type.append(ws.cell_value(1,3))
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
    file_object.close()

    #reading the rest of the files
    for i in range(2, ws.nrows):
        name = ws.cell_value(i,0)
        sample_type.append(ws.cell_value(i,3))
        file_object = open(datapath+"data"+tissue_id+"/"+name, "r")
        lines = file_object.readlines()
        if(ws.cell_value(i,3) == 'Solid Tissue Normal'):
            normal.append(i - 1)
        else:
            tumor.append(i - 1)
        data.append([])
        for line in lines:
            data[-1].append(line.split()[1])
        file_object.close()

    data = np.array(data, dtype="f") + 0.1
    ref = gmean(data[normal])
    data = np.log2(data/ref)

    return([data, normal, tumor])


# Reading and processing the data ------------------------------------------

print("Reading files from databases:")
data, normal, tumor = read_data(datapath,tissue_id)
print("Data successfully loaded!")
print("normal =", len(normal))
print("tumor =", len(tumor))
print("total =", len(data))

# reducing the data by Ngen genes for Lite version
if(lite_version and data.shape[1] > Ngen):
    print("Lite version: reducing the number of genes by", Ngen)
    data = data[:, :Ngen]
elif(len(data) < Ngen):
    lite_version = False
    print("Warning: Imposible reduction, number of genes larger than databases.")
    print("The Full version will be performe instead.")

print("Computing PCA components...")
if(lite_version):
    eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PC_decomp(data,normal)
else:
    eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PC_decomp(data,normal,True)

principal = eigenvectors[:, eigenvalues.argmax()]
index = np.argpartition(-np.abs(principal), 20)[:20]
components = principal[index]
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
    np.savetxt('ind-normal'+tissue_id+'_dat.dat', normal, fmt='%i')
    np.savetxt('ind-tumor'+tissue_id+'_dat.dat', tumor, fmt='%i')

    plt.scatter(projection[0,normal], -projection[1,normal], c='b', s=15,label="Normal")
    plt.scatter(projection[0,tumor], -projection[1,tumor], c='r', s=15,label="Tumor")
    plt.grid(True)
    plt.xlabel('PC1')
    plt.ylabel('-PC2')
    plt.legend()
    plt.savefig(tissue_id+'_PC2_Vs_PC1_fig.png')
else:
    np.savetxt(outputpath+'pc'+tissue_id+'.xls', projection, fmt='%f')

print("Done!")
