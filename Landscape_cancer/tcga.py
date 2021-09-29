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
from Manifolds_normal_cancer import PCA_core_mem as pca_mem
from os.path import *


# General variables ------------------------------------------------------

lite_version = True
Ngen = 1000 #60483
NgenS = 20 #must be smaller or equal to Npc
Npc = 30; NpcL = 100 #the maximum possible number of PCs is the number of genes
tissue_id = "KIRC"
#possible targets:
#tissue_id=["PRAD", "KIRC", "LUSC", "LUAD", "UCEC", "KIRP", "BLCA", "COAD", 
#                   "ESCA", "LIHC", "STAD", "THCA", "BRCA", "HNSC", "READ"]

datapath = "../databases_external/TCGA/"
outputpath = "../databases_generated/TCGA_pca/"  # for Full databases


# Functions ------------------------------------------------------

def read_GEdata(datapath,tissue_id):
    wb = xl.open_workbook(datapath+"sample"+tissue_id+".xls")
    ws = wb.sheet_by_index(0)
    print(ws.nrows,"individual files were identified.")

    sample_type = []
    normal = []
    tumor = []
    GEdata = []
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
    GEdata.append([])
    for line in lines:
        GEdata[-1].append(line.split()[1])
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
        GEdata.append([])
        for line in lines:
            GEdata[-1].append(line.split()[1])
        file_object.close()

    GEdata = np.array(GEdata, dtype="f") + 0.1
    ref = gmean(GEdata[normal])
    refT = gmean(GEdata[tumor,:])/ref
    GEdata = np.log2(GEdata/ref)

    refTu = []; refTo = []
    for i in range(len(refT)):
        if refT[i]<1.:
            refTu.append(refT[i])
        else:
            refTo.append(refT[i])
    refTu.sort()
    refTo.sort()

    return([GEdata, refTu, refTo, normal, tumor, np.array(genes_id)])


# Reading and processing the data ------------------------------------------

print("Reading files from "+tissue_id+" databases:")
if(not isfile(datapath+"sample"+tissue_id+".xls")):
    print(tissue_id, "databases not found, please read", "'"+ datapath + "info.txt'")
    sys.exit()
data, refTu, refTo, normal, tumor, genes_id = read_GEdata(datapath,tissue_id)
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
    #eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca_mem.PCAL(data,NpcL, matrix_fname='tri_tcga_8000_dat.npy', buffer_size=8000)

    radius_pc1_normal, center_pc1_normal = pca.region(projection[0,normal])
    radius_pc1_tumor, center_pc1_tumor = pca.region(projection[0,tumor])
    radius_pc2_normal, center_pc2_normal = pca.region(projection[1,normal])
    radius_pc2_tumor, center_pc2_tumor = pca.region(projection[1,tumor])
    theta = np.linspace(0, 2*np.pi, 100)

else:
    eigenvalues,eigenvectors,eigenvalues_normalized,projection = pca.PCA_core(data,True,Npc)

index = np.argsort(-np.abs(eigenvectors[:, 0]))[:Npc]
PC1 = eigenvectors[index, 0]
print("Done!")


# Output data and figures --------------------------------------------------

print("Exporting data and plots...")
if(lite_version):
    #np.savetxt('pc'+tissue_id+'_dat.dat', projection.T, fmt='%f')
    #np.savetxt('evec'+tissue_id+'_dat.dat', eigenvectors.T, fmt='%f')
    np.savetxt(tissue_id+'eigenValN_dat.dat', eigenvalues_normalized, fmt='%f', header='PCA normalized eigenvalues for '+tissue_id)
    #np.savetxt('evec-t'+tissue_id+'_dat.dat', eigenvectors, fmt='%f')
    #np.savetxt('ind20'+tissue_id+'_dat.dat', index, fmt='%i')
    np.savetxt(tissue_id+'geneID-PC1_dat.dat', np.transpose([index,genes_id[index],PC1]), 
    delimiter="\t\t", fmt="%s", header='index\tgene ID\t\t\tPC1')

   #Fig.1
    fig1, ax1 = plt.subplots(1,2)
    f=0.95  
    fig1.set_size_inches(f*2*4, f*1*3)
    plt.rcParams['lines.linewidth'] = 1.0
    plt.rcParams['axes.labelsize'] = 15
    plt.rcParams['xtick.labelsize'] = 15
    plt.rcParams['ytick.labelsize'] = 15
    textF='large'
    colors = ['royalblue','r','tab:orange','forestgreen','k']
    al=0.7
    size=18
    xmin=-6;xmax=35
    ymin=-21;ymax=45
    
    ax1[0].plot([xmin,xmax],[0, 0],c='gray',linestyle='--')
    ax1[0].plot([0, 0],[ymin,ymax],c=colors[0])
    ax1[0].plot([center_pc1_tumor, center_pc1_tumor],[ymin,ymax],c=colors[2])
    ax1[0].scatter(projection[0,normal], -projection[1,normal],s=size,label=tissue_id+', normal', marker="o",c='b',alpha=al)
    ax1[0].scatter(projection[0,tumor], -projection[1,tumor],s=size,label=tissue_id+', tumor', marker="s",c='r',alpha=al)
    ax1[0].plot(center_pc1_normal+radius_pc1_normal*np.cos(theta),center_pc2_normal+radius_pc2_normal*np.sin(theta),c=colors[0])
    ax1[0].plot(center_pc1_tumor+radius_pc1_tumor*np.cos(theta),center_pc2_tumor+radius_pc2_tumor*np.sin(theta),c=colors[2])
    ax1[1].plot(range(1,NgenS+1),100*np.cumsum(eigenvalues_normalized[:NgenS]), marker="o",c='b',markersize=5)
    ax1[1].text(3, 25, tissue_id, size=textF)

    ax1[0].set_xlabel('PC1'), ax1[0].set_ylabel('-PC2'), ax1[0].legend(edgecolor='k', borderaxespad=0.2)
    ax1[0].locator_params(axis='x',nbins=6),ax1[0].locator_params(axis='y',nbins=6), ax1[0].set_xlim(xmin,xmax), ax1[0].set_ylim(ymin,ymax)
    ax1[1].set_xlabel('Number of components'), ax1[1].set_ylabel('Cumulative variance (%)')
    ax1[1].set_xlim(0.,NgenS+0.2), ax1[1].set_ylim(0,102), ax1[1].locator_params(axis='x',nbins=5)  
    fig1.tight_layout()

   #Fig.2
    fig2, ax2 = plt.subplots(1,3)  
    f=0.99 
    fig2.set_size_inches(f*3*4, f*1*3) 
    plt.rcParams['lines.linewidth'] = 1.5
    plt.rcParams['axes.labelsize'] = 5
    plt.rcParams['xtick.labelsize'] = 5
    plt.rcParams['ytick.labelsize'] = 5
    textF='large'
    colors = ['royalblue','r','tab:orange','forestgreen','k']
    al=0.7
    size=15

    ax2[0].plot([0.,(Ngen+0.5)/1000.],[0, 0],c='gray',linestyle='--')
    for i in index:
        ax2[0].plot([(i+1)/1000.,(i+1)/1000.],[0, eigenvectors[i,0]], c='b')
    Ng_maxU=1000;Ng_maxO=1000
    ax2[1].loglog(refTu[:Ng_maxU],range(1,Ng_maxU+1), mec='r', markersize=5, marker='o',mfc='w', ls='')
    ax2[2].loglog(refTo[-Ng_maxO:],range(Ng_maxO,0,-1), mec='r', markersize=5, marker='o',mfc='w', ls='')

    ax2[0].set_xlabel('Gene number ('+r'$\times$1000)'), ax2[0].set_ylabel('Amplitude') 
    ax2[1].set_xlabel('Gene expression'), ax2[1].set_ylabel('Number of genes'), ax2[2].set_xlabel('Gene expression')
    ax2[0].set_xlim(0.,(Ngen+0.5)/1000.)#,ax2[1].set_xlim(0.5*10**-4,1.5)
    ax2[0].locator_params(axis='x',nbins=6),ax2[1].set_xticks([10**-4, 10**-3, 10**-2,10**-1,1]),ax2[2].set_xticks([1, 10, 10**2,10**3])
    ax2[1].grid(which='major', ls='-'); ax2[2].grid(which='major', ls='-')
    ax2[1].grid(which='minor', ls='--'); ax2[2].grid(which='minor', ls='--')
    fig2.tight_layout() 
      
    fig1.savefig(tissue_id+'_PC1-2_fig.pdf')
    fig2.savefig(tissue_id + "_gene-dist_fig.pdf")
else:
    np.savetxt(outputpath+'pc'+tissue_id+'.xls', projection.T)

print("Done!")
