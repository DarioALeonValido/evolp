"""
This script file processes PCA-TCGA data of a set of 
tissues. The data sources are mentioned in the 
imported files, while the procedures and main results 
are described in paper arxiv:2003.07828v1

For more details on this file see author(s):
DALV
"""

import numpy as np
import xlrd as xl
import matplotlib as mpl
import matplotlib.pyplot as plt


# General variables -----------------------------------------------

NpcP = 2
NpcV = 18
initpc = 1
tissue_id = 'KIRP'
tissues_examples = ['KIRP','LIHC','LUSC']
tissues_PCA = ['BLCA','BRCA','COAD','ESCA','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','READ','STAD','THCA','UCEC']

sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'


# Functions ------------------------------------------------------

def read_PC(sample_path,PC_path,tissue_id,initpc,Npc):
    ind_normal = []
    ind_tumor = []
    pc_data = []
    sample_wb = xl.open_workbook(sample_path+"sample"+tissue_id+".xls")
    sample_ws = sample_wb.sheet_by_index(0)
    pc_wb = xl.open_workbook(PC_path+"pc"+tissue_id+".xls")
    pc_ws = pc_wb.sheet_by_index(0)
    print('Number of',tissue_id+' cases:', sample_ws.nrows)
    for i in range(sample_ws.nrows-1):
        if(sample_ws.cell_value(i+1,3) == 'Solid Tissue Normal'):
            ind_normal.append(i)
        else:
            ind_tumor.append(i)
        pc_data.append([])
        for j in range(initpc-1,initpc-1+Npc):    
            pc_data[i].append(pc_ws.cell_value(i,j))

    pc_data = np.array(pc_data, dtype="f")

    return pc_data, ind_normal, ind_tumor


def fitness_dist(pc1,ind_n,ind_t,fitN,fitT,binsN,binsT):

    histN,bin_edgesN = np.histogram(pc1[ind_n],bins=binsN)
    histT,bin_edgesT = np.histogram(pc1[ind_t],bins=binsT)
    pc = np.concatenate((bin_edgesN[0:-1],bin_edgesT[0:-1]), axis=None)
    mfitness = np.concatenate((-fitN*histN/max(histN),-fitT*histT/max(histT)))

    return pc,mfitness


def compute_GEdistances(sample_path,PC_path,tissues,ipc):
    Xt = []
    Rn = []
    Rt = []
    #D = []
    for t_id in tissues:
        pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,t_id,ipc,1)
        Xn = np.mean(pc_data[ind_normal])
        Rn.append(np.std(pc_data[ind_normal],ddof=1))
        Xt.append(abs(np.mean(pc_data[ind_tumor]) - Xn))
        Rt.append(np.std(pc_data[ind_tumor]-Xn,ddof=1))  
    Xt = np.array(Xt, dtype="f")      
    Rn = np.array(Rn, dtype="f")   
    Rt = np.array(Rt, dtype="f")   

    return Xt,Rn,Rt#,D


def compute_cum_variance(sample_path,PC_path,tissue,Npc): # this function is useless
    CV  = []
    CVn = []
    CVt = []
    pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissue,1,1)
    CVn.append(np.std(pc_data[ind_normal],ddof=1))
    CVt.append(np.std(pc_data[ind_tumor],ddof=1))
    CV.append(np.std(pc_data,ddof=1))
    for i in range(1,Npc):
        pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissue,i,1)
        CVn.append(CVn[-1] + np.std(pc_data[ind_normal],ddof=1))
        CVt.append(CVt[-1] + np.std(pc_data[ind_tumor],ddof=1))
        CV.append(CV[-1] + np.std(pc_data,ddof=1)) 
    CV  = np.array(CV, dtype="f")      
    CVn = np.array(CVn, dtype="f")   
    CVt = np.array(CVt, dtype="f")   

    return CV,CVn,CVt


# Reading and processing all the data ---------------------------------------

print("\nLoading Principal Componets for:") #working with PCA data
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_PCA,initpc) # some plots do not work with i = 2, 3, ...
R=Xt-(Rt+Rn)


# Exporting readable data ------------------------------------------------------

with open('geData_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of GE for each tissue.\n\n')
    savefile.write('Tissue\tXt\t\tRn\t\tRt\t\tR\n')
    for i in range(len(tissues_PCA)):
        savefile.write(' '+tissues_PCA[i]+'\t')
        savefile.write('%f\t' % Xt[i])
        savefile.write('%f\t' % Rn[i])
        savefile.write('%f\t' % Rt[i])
        savefile.write('%f\n' % R[i])
    savefile.write('\n')
    savefile.close()


# Plots ---------------------------------------------------------------

fig1, axs1 = plt.subplots(3, 2, sharex=True)
f=0.95
fig1.set_size_inches(f*3*3, f*2*4)
fig1.subplots_adjust(wspace=0.001, hspace=0.0001)
#plt.rcParams['font.size'] = '12'
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.labelsize'] = 15
plt.rcParams['xtick.labelsize'] = 15
plt.rcParams['ytick.labelsize'] = 15
textF='large'

colors = ['royalblue','r','tab:orange','forestgreen','k']
al=0.7
size=18
  
# plotting PC1 vs. PC2 ----------------------------
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_examples,initpc) # do not currently work with i = 2, 3, ...
xmin=-50;xmax=285
ymin=-230;ymax=180
label=''
for t in range(len(tissues_examples)):
    pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissues_examples[t],initpc,NpcP)
    axs1[t,0].plot([xmin,xmax],[0, 0],c='gray',linestyle='--')
    axs1[t,0].plot([0, 0],[ymin,ymax],c='b')
    axs1[t,0].plot([Xt[t], Xt[t]],[ymin,ymax],c='r')
    #CV,CVn,CVt = compute_cum_variance(sample_path,PC_path,tissues_examples[t],NpcV)
    #print(CVn)
    if (tissues_examples[t]=='KIRP'):                                                                     #, facecolors='none',edgecolors='b'
        axs1[t,0].scatter(-pc_data[ind_normal,0],pc_data[ind_normal,1],label=tissues_examples[t]+', normal',s=size,c='b',alpha=al)
        axs1[t,0].scatter(-pc_data[ind_tumor,0],pc_data[ind_tumor,1],label=tissues_examples[t]+', tumor',s=size, marker="s",c='r',alpha=al)
        pc,mfitness = fitness_dist(-pc_data[:,0],ind_normal,ind_tumor,1.,1.5,5,16)
        axs1[t,1].plot(pc,mfitness,label='-fitness',c='black')
    else:
        axs1[t,0].scatter(pc_data[ind_normal,0],-pc_data[ind_normal,1],label=tissues_examples[t]+', normal',s=size,c='b',alpha=al)
        axs1[t,0].scatter(pc_data[ind_tumor,0],-pc_data[ind_tumor,1],label=tissues_examples[t]+', tumor',s=size, marker="s",c='r',alpha=al)
        pc,mfitness = fitness_dist(pc_data[:,0],ind_normal,ind_tumor,1.,1.5,5,16)
        axs1[t,1].plot(pc,mfitness,label='-fitness',c='black')

    axs1[t,0].set_xlim(xmin,xmax)
    axs1[t,0].set_ylim(ymin,ymax)
    axs1[t,0].set_ylabel("PC2")
    #axs1[t,0].text(-40, -180, tissues_examples[t], size='large')
    axs1[t,0].legend(bbox_to_anchor=(0., 0.),edgecolor='k', loc='lower left', borderaxespad=0.2)
    axs1[t,0].locator_params(axis='x', nbins=7)
    axs1[t,0].locator_params(axis='y', nbins=9)


    axs1[t,1].set_ylabel("-Fitness")
    axs1[t,1].text(0, -1.2, tissues_examples[t])
    axs1[t,1].locator_params(axis='x', nbins=7)

    if t != len(tissues_examples)-1:
        label = label + tissues_examples[t] + '-'
    else:
        label = label + tissues_examples[t]   
    
axs1[t,0].set_xlabel("PC1")
axs1[t,1].set_xlabel("PC1")
plt.tight_layout()
plt.savefig('ge'+label+'_pc1-2_fig.pdf')
plt.clf()


print("Done!\n")
print("Check out the outputs!\n")
