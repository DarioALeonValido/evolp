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
from matplotlib.ticker import AutoMinorLocator
import seaborn as sns
import math
#import plotly.graph_objects as go


# General variables -----------------------------------------------

Ngen = 60483
NpcP = 2
NpcV = 18
initpc = 1
tissue_id = 'KIRC'
tissues_examples = ['KIRP','LIHC','LUSC']
genes_examples=['UMOD','CYPIA2','SFTPC']
tissues_PCA = ['BLCA','BRCA','COAD','ESCA','HNSC','KIRC','KIRP','LIHC','LUAD','LUSC','PRAD','READ','STAD','THCA','UCEC']

sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'
GE_path = '../databases_generated/TCGA_exp/'
td_path = './data/'
commG_file = 'common-genes.dat'


# Functions ------------------------------------------------------

def read_eigenval(tissue_id,PC_path):
    file = open(PC_path+'eigenvalues'+tissue_id+'.dat', 'r')
    lines = file.readlines()

    eigenval = []     
    for line in lines[1:]:
        eigenval.append(line.split()[0])
    eigenval = np.array(eigenval, np.float)
  
    return eigenval 

def read_genesA(tissue_id,PC_path):
    file = open(PC_path+'genes'+tissue_id+'.dat', 'r')
    lines = file.readlines()

    indexes = []; amplitudes = []     
    for line in lines[1:]:
        indexes.append(line.split()[0])
        amplitudes.append(line.split()[-1])
    indexes = np.array(indexes, np.int)
    amplitudes = np.array(amplitudes, np.float)
  
    return indexes,amplitudes 

def read_EXPdist(tissue_id,GE_path,over=True):
    if over:
        exp = 'over'
    else:
        exp = 'under'
    file = open(GE_path+'dist-'+exp+tissue_id+'.dat', 'r')
    lines = file.readlines()

    exp = []     
    for line in lines[1:]:
        exp.append(line.split()[0])
    if over:
        exp = np.array(exp, np.float)
    else:
        exp = np.array(exp, np.float)
  
    return exp 

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

def find_stage(ID, clinical_ws):
    for i in range(1, clinical_ws.nrows):
        if clinical_ws.cell_value(i,0) == ID:
            return clinical_ws.cell_value(i,4)

def same_stage(stage, s):
    same_stage = False
    if (stage=='Stage '+s or stage=='Stage '+s+'A' or stage=='Stage '+s+'B' or stage=='Stage '+s+'C'):
        same_stage = True

    return same_stage

def read_PC_clinical(sample_path,PC_path,tissue_id,initpc,Npc):
    ind_normal = []
    ind_tumor = []
    pc_data = []
    sample_wb = xl.open_workbook(sample_path+"sample"+tissue_id+".xls")
    sample_ws = sample_wb.sheet_by_index(0) 
    pc_wb = xl.open_workbook(PC_path+"pc"+tissue_id+".xls")
    pc_ws = pc_wb.sheet_by_index(0)

    ind_stageI = []
    ind_stageII = []
    ind_stageIII = []
    ind_stageIV = []
    clinical_wb = xl.open_workbook(sample_path+"clinical"+tissue_id+".xls")
    clinical_ws = clinical_wb.sheet_by_index(0)

    #print('Number of',tissue_id+' cases:', sample_ws.nrows)
    for i in range(sample_ws.nrows-1):
        if(sample_ws.cell_value(i+1,3) == 'Solid Tissue Normal'):
            ind_normal.append(i)
        else:
            ind_tumor.append(i)
            #if find_stage(sample_ws.cell_value(i+1,1), clinical_ws) == 'Stage I': ind_stageI.append(i)
            if same_stage(find_stage(sample_ws.cell_value(i+1,1), clinical_ws),'I'): ind_stageI.append(i)
            elif same_stage(find_stage(sample_ws.cell_value(i+1,1), clinical_ws),'II'): ind_stageII.append(i)
            elif same_stage(find_stage(sample_ws.cell_value(i+1,1), clinical_ws),'III'): ind_stageIII.append(i)
            elif same_stage(find_stage(sample_ws.cell_value(i+1,1), clinical_ws),'IV'): ind_stageIV.append(i)            
        pc_data.append([])
        for j in range(initpc-1,initpc-1+Npc):    
            pc_data[i].append(pc_ws.cell_value(i,j))
    pc_data = np.array(pc_data, dtype="f")

    return pc_data, ind_normal, ind_tumor, ind_stageI, ind_stageII, ind_stageIII, ind_stageIV

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

def pc_sign(pc_array):
    sign = np.sign(np.mean(pc_array))
    #sign = np.sign(pc_array[0])
    #maxvalue = abs(pc_array[0])
    #for i in range(1,len(pc_array)):
    #    if abs(pc_array[i])>maxvalue:
    #        sign = np.sign(pc_array[i])
    #        maxvalue = abs(pc_array[i])

    return sign#, sign*pc_array


def read_comm_genes(gfile):
    file = open(gfile, 'r')
    lines = file.readlines()

    Ngen = int(lines[0].split()[-1])
    tissues = []; corr_mat = []  
    tissues.append(lines[2].split()[0]) 
    corr_mat.append([Ngen])  
    for i in range(3,len(lines)):
        tissues.append(lines[i].split()[0])
        corr_mat.append([])
        for j in range(1,i-1):
            corr_mat[-1].append(int(lines[i].split()[j]))
            corr_mat[i-2-j].append(int(lines[i].split()[i-1-j]))
        corr_mat[i-2].append(Ngen)
    #corr_mat = np.array(corr_mat, np.int)  

    return tissues,corr_mat 

def read_tot_density(tissue_id,td_path):
    file = open(td_path+tissue_id+'_total-density.dat', 'r')
    lines = file.readlines()

    axis_pc1 = []
    axis_pc2 = []
    tot_density = []    
    axis_pc1.append(float(lines[0].split()[0])) 
    axis_pc2.append(float(lines[0].split()[1]))
    tot_density.append([]) 
    tot_density[-1].append(float(lines[0].split()[2]))
    inter = 1
    while float(lines[inter].split()[0]) == axis_pc1[-1]:
        axis_pc2.append(float(lines[inter].split()[1])) 
        tot_density[-1].append(float(lines[inter].split()[2]))
        inter=inter+1
    for i in range(inter,len(lines)):
        if (i)%inter==0:
            axis_pc1.append(float(lines[i].split()[0]))
            tot_density.append([])
        tot_density[-1].append(float(lines[i].split()[2]))
    axis_pc1 = np.array(axis_pc1, np.float)
    axis_pc2 = np.array(axis_pc2, np.float)
    tot_density = np.array(tot_density)

    return axis_pc1, axis_pc2, tot_density.T 

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


# Reading and processing all the data ---------------------------------------

print("\nLoading Principal Components for:") #working with PCA data
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_PCA,initpc) # some plots do not work with i = 2, 3, ...
R=Xt-(Rt+Rn)

comm_tiss,corr_mat = read_comm_genes(PC_path+commG_file)

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


# Settings for fig1 ---------------------------------------------------------------
plt.rcParams['font.size'] = 11.5
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['lines.markersize'] = 5.2

fig1, axs1 = plt.subplots(3, 2, sharex='col')
f=0.95
fig1.set_size_inches(f*3*3, f*2*4.5)
fig1.subplots_adjust(wspace=0.001, hspace=0.0001)
colors = ['royalblue','r','tab:orange','forestgreen','k']
al=0.7

# plotting PC1 vs. PC2 ----------------------------
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_examples,initpc) # do not currently work with i = 2, 3, ...
xmin=-50;xmax=285
ymin=-230;ymax=180
sign1=[];sign2=[]
label=''
for t in range(len(tissues_examples)):
    pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissues_examples[t],initpc,NpcP)
    axs1[t,0].plot([xmin,xmax],[0, 0],c='gray',linestyle='--')
    axs1[t,0].plot([0, 0],[ymin,ymax],c='b')
    axs1[t,0].plot([Xt[t], Xt[t]],[ymin,ymax],c='r')
    eigenval = read_eigenval(tissues_examples[t],PC_path)

    sign1.append(pc_sign(pc_data[ind_tumor,0]))
    sign2.append(pc_sign(pc_data[ind_tumor,1]))                                                  #, facecolors='none',edgecolors='b'
    axs1[t,0].scatter(sign1[-1]*pc_data[ind_normal,0],sign2[-1]*pc_data[ind_normal,1],label=tissues_examples[t]+', normal',c='b',alpha=al)
    axs1[t,0].scatter(sign1[-1]*pc_data[ind_tumor,0],sign2[-1]*pc_data[ind_tumor,1], 
    label=tissues_examples[t]+', tumor', marker="s",c='r',alpha=al)
    #pc,mfitness = fitness_dist(-pc_data[:,0],ind_normal,ind_tumor,1.,1.5,5,16)
    #axs1[t,1].plot(pc,mfitness,label='-fitness',c='black')
    axs1[t,1].plot(range(1,NpcV+1),100*np.cumsum(eigenval[:NpcV]), marker="o",c='b')

    axs1[t,0].set_xlim(xmin,xmax)
    axs1[t,0].set_ylim(ymin,ymax)
    axs1[t,0].set_xticks([-50,0,50,100,150,200,250])
    axs1[t,0].set_yticks([-200,-100,0,100])
    axs1[t,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs1[t,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    axs1[t,0].set_ylabel("PC2")
    #axs1[t,0].text(-40, -180, tissues_examples[t])
    axs1[t,0].legend(bbox_to_anchor=(0., 0.),edgecolor='k', loc='lower left', borderaxespad=0.2)
    #axs1[t,0].locator_params(axis='x', nbins=7)
    #axs1[t,0].locator_params(axis='y', nbins=9)


    axs1[t,1].set_ylabel("Cumulative variance (%)")
    axs1[t,1].text(5, 20, tissues_examples[t])
    #axs1[t,1].locator_params(axis='x', nbins=7)
    axs1[t,1].set_xticks([0,2,4,6,8,10,12,14,16,18])
    axs1[t,1].set_yticks([0,20,40,60,80,100])
    axs1[t,1].xaxis.set_minor_locator(AutoMinorLocator(2)),axs1[t,1].yaxis.set_minor_locator(AutoMinorLocator(2))

    if t != len(tissues_examples)-1:
        label = label + tissues_examples[t] + '-'
    else:
        label = label + tissues_examples[t]
      
axs1[t,0].set_xlabel("PC1")
axs1[t,1].set_xlabel("Number of components")
plt.tight_layout()
plt.savefig(label+'_pc1-2_fig.pdf')
#plt.clf()


# Settings for fig2 ---------------------------------------------------------------
#mpl.style.use
plt.rcParams['font.size'] = 12.
plt.rcParams['lines.linewidth'] = 1.5
plt.rcParams['axes.labelsize'] = 12.
plt.rcParams['xtick.labelsize'] = 11.9
plt.rcParams['ytick.labelsize'] = 11.9
plt.rcParams['lines.markersize'] = 5.0

fig2, axs2 = plt.subplots(3, 3)
#fig2.subplots_adjust(wspace=0.0, hspace=1.8)
#fig2.subplots_adjust(left=0.0, bottom=0.0, right=0.9, top=0.9, wspace=0., hspace=0.1)
f=0.99
fig2.set_size_inches(f*3*4, f*3*3)
#fig2.subplots_adjust(wspace=0.001, hspace=0.0001)
colors = ['royalblue','r','tab:orange','forestgreen','k']
al=0.7
xmin=0;xmax=63
ymin=-0.08;ymax=0.05
  
# plotting gene distribution ----------------------------
#label=''
for t in range(len(tissues_examples)):
    pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissues_examples[t],initpc,NpcP)
    #axs2[t,0].plot([0.,Ngen/1000.],[0, 0],c='gray',linestyle='--')
    axs2[t,0].plot([xmin,xmax],[0, 0],c='gray',linestyle='--')
    indexes,amplitudes = read_genesA(tissues_examples[t],PC_path)
    ip=0;peak=abs(amplitudes[0]) 
    for i in range(len(indexes)):
        axs2[t,0].plot([(indexes[i]+1)/1000.,(indexes[i]+1)/1000.],[0, sign1[t]*amplitudes[i]], c='b')
        if abs(amplitudes[i])>peak:
            ip=i;peak=abs(amplitudes[i])
    axs2[t,0].plot([(indexes[ip]+1)/1000.,(indexes[ip]+1)/1000.],[0, sign1[t]*amplitudes[ip]], c='r')
    #axs2[t,0].scatter([(indexes[ip]+1)/1000.],[pc_sign(pc_data[ind_normal,0])[0]*amplitudes[ip]], c='r', label=genes_examples[t])
    axs2[t,0].text((indexes[ip]+1)/1000.+1.5,sign1[t]*amplitudes[ip]-0.0025,genes_examples[t],ha='left',c='r')

    expO = read_EXPdist(tissues_examples[t],GE_path,True)
    expU = read_EXPdist(tissues_examples[t],GE_path,False)
    axs2[t,1].loglog(expU,range(1,len(expU)+1), mec='r', marker='o',mfc='w', ls='',label=tissues_examples[t]+', under-exp')
    axs2[t,2].loglog(expO,range(len(expO),0,-1), mec='r', marker='o',mfc='w', ls='',label=tissues_examples[t]+', over-exp')


    axs2[t,0].set_xlim(xmin,xmax)
    #axs1[t,0].set_ylim(ymin,ymax)
    axs2[t,0].set_ylabel('Amplitude, '+tissues_examples[t]), axs2[t,1].set_ylabel('Number of genes')
    #axs2[t,0].set_yticks([-0.08,-0.06,-0.04,-0.02,0.,0.02,0.04,0.06])
    axs2[t,0].locator_params(axis='y',nbins=6)
    axs2[t,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs2[t,0].yaxis.set_minor_locator(AutoMinorLocator(2))
    axs2[t,1].set_xticks([10**-4, 10**-3, 10**-2,10**-1,1]),axs2[t,2].set_xticks([1, 10, 10**2,10**3,10**4])
    axs2[t,2].tick_params(labelleft=False)
    #axs2[t,1].text(2*10**-4,6*10**2, tissues_examples[t]+', under-exp',ha='left')
    #axs2[t,2].text(3*10**1,6*10**2, tissues_examples[t]+', over-exp',ha='left')
    axs2[t,1].grid(which='major', axis='both', ls='-'),axs2[t,1].grid(which='minor', axis='both', ls='--')
    axs2[t,2].grid(which='major', axis='both', ls='-'),axs2[t,2].grid(which='minor', axis='both', ls='--')
    #axs2[t,0].legend(bbox_to_anchor=(0., 0.),edgecolor='k', loc='lower left', borderaxespad=0.1)
    axs2[t,1].legend(bbox_to_anchor=(0., 1.),edgecolor='k', loc='upper left', borderaxespad=0.4)
    axs2[t,2].legend(bbox_to_anchor=(1., 1.),edgecolor='k', loc='upper right', borderaxespad=0.4)

#    if ti != len(tissues_examples)-1:
#        label = label + tissues_examples[ti] + '-'
#    else:
#        label = label + tissues_examples[ti]   

axs2[0,0].tick_params(labelbottom=False),axs2[0,1].tick_params(labelbottom=False),axs2[0,2].tick_params(labelbottom=False)
axs2[1,0].tick_params(labelbottom=False),axs2[1,1].tick_params(labelbottom=False),axs2[1,2].tick_params(labelbottom=False)
axs2[t,0].set_xlabel('Gene number ('+r'$\times$1000)')
axs2[t,1].set_xlabel('Differential expression'), axs2[t,2].set_xlabel('Differential expression')
plt.tight_layout()
plt.savefig(label+'_genes_fig.pdf')
#plt.clf()

# Settings for fig table 2 ---------------------------------------------------------------
mpl.style.use
plt.rcParams['xtick.labelsize'] = 5.5
plt.rcParams['ytick.labelsize'] = 5.5
heat_scale = 2000
for i in range(len(corr_mat)): corr_mat[i][i]=heat_scale
mask = np.zeros_like(corr_mat, dtype=np.bool)
mask[np.triu_indices_from(mask)]= True
sns.set_style({'xtick.bottom': True}, {'ytick.left': True})

# plotting fig table 2 ----------------------------
fig3, ax = plt.subplots()
heatmap = sns.heatmap(corr_mat, fmt='d', mask = mask, square = True, linewidths = .2, cmap = 'YlOrBr',
                      cbar_kws = {#'orientation': 'horizontal', #'shrink': .4,
                                'ticks' : [0, heat_scale/4, heat_scale/2, heat_scale*3/4, heat_scale]},
                      vmin = 0, vmax = heat_scale, annot = True, annot_kws = {'size': 5.5})

bottom, top = ax.get_ylim()
ax.set_ylim(bottom + 0.5, top - 0.5)
ax.set_yticklabels(comm_tiss, rotation = 0)
ax.set_xticklabels(comm_tiss)
heatmap.get_figure().savefig('common-genes_fig.pdf', bbox_inches='tight')
#plt.clf()

# Settings for fig 4 ---------------------------------------------------------------
mpl.style.use
plt.rcParams['lines.linewidth'] = 1.0
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11
plt.rcParams['legend.fontsize'] = 11
plt.rcParams['lines.markersize'] = 4.

# plotting fig 4 ----------------------------
fig4, axs4 = plt.subplots(2,2,sharex='all', sharey='all')
f=0.95  
fig4.set_size_inches(f*2*4, f*2*3.8)
al=0.7

axis_pc1, axis_pc2, tot_density = read_tot_density(tissue_id,td_path)
grid_pc1, grid_pc2 = np.meshgrid(axis_pc1, axis_pc2)

pc_data,ind_normal,ind_tumor,ind_stageI,ind_stageII,ind_stageIII,ind_stageIV = read_PC_clinical(sample_path,PC_path,tissue_id,initpc,NpcP)
sign1=pc_sign(pc_data[ind_tumor,0]);sign2=-pc_sign(pc_data[ind_tumor,1])

colormap='pink_r'#'autumn_r'#'pink_r'#'gist_heat_r'#'OrRd'#'RdYlGn_r'#'RdYlBu_r'#'Oranges'#'YlOrBr'#'Greys'
npoints=20
axs4[0,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4[0,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4[1,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4[1,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)

axs4[0,0].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4[0,0].scatter(sign1*pc_data[ind_stageI,0],sign2*pc_data[ind_stageI,1],label=tissue_id+', stage I', marker="s",c='r',alpha=al)
axs4[0,1].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4[0,1].scatter(sign1*pc_data[ind_stageII,0],sign2*pc_data[ind_stageII,1],label=tissue_id+', stage II', marker="s",c='r',alpha=al)
axs4[1,0].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4[1,0].scatter(sign1*pc_data[ind_stageIII,0],sign2*pc_data[ind_stageIII,1],label=tissue_id+', stage III', marker="s",c='r',alpha=al)
axs4[1,1].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4[1,1].scatter(sign1*pc_data[ind_stageIV,0],sign2*pc_data[ind_stageIV,1],label=tissue_id+', stage IV', marker="s",c='r',alpha=al)

axs4[1,0].set_xlabel('PC1'), axs4[1,1].set_xlabel('PC1')
axs4[0,0].set_ylabel('PC2'), axs4[1,0].set_ylabel('PC2')
axs4[0,0].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4[0,1].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4[1,0].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4[1,1].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4[0,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4[0,0].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4[0,1].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4[0,1].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4[1,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4[1,0].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4[1,1].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4[1,1].yaxis.set_minor_locator(AutoMinorLocator(2))
fig4.tight_layout()
fig4.savefig(tissue_id + "_stages_fig.pdf")

# plotting fig 4 with generated data ----------------------------
tissue_id = 'LIHC'
fig4g, axs4g = plt.subplots(2,2,sharex='all', sharey='all')
f=0.95  
fig4g.set_size_inches(f*2*4, f*2*3.8)
al=0.7

dev=30.

#KIRC
dev=[40.,50.]
ranges=[[-50.,300.],[-140.,430.]]
npoints=[36,58]

#KIRP&LUSC
dev=[40.,50.]
ranges=[[-50.,300.],[-180.,250.]]
npoints=[36,44]

#LIHC
dev=[40.,50.]
ranges=[[-80.,300.],[-180.,250.]]
npoints=[39,44]

pc_data,ind_normal,ind_tumor,ind_stageI,ind_stageII,ind_stageIII,ind_stageIV = read_PC_clinical(sample_path,PC_path,tissue_id,initpc,NpcP)
sign1=pc_sign(pc_data[ind_tumor,0]);sign2=-pc_sign(pc_data[ind_tumor,1])

axis_pc1, axis_pc2, tot_density = compute_density(np.transpose([sign1*pc_data[:,0],sign2*pc_data[:,1]]),dev,ranges,npoints)
grid_pc1, grid_pc2 = np.meshgrid(axis_pc1, axis_pc2)

dev=[25.,40.]
pc1_normal, pc2_normal, density_normal = compute_density(np.transpose([sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1]]),dev,ranges,npoints)
dev=[25.,40.]
pc1_tumor, pc2_tumor, density_tumor = compute_density(np.transpose([sign1*pc_data[ind_tumor,0],sign2*pc_data[ind_tumor,1]]),dev,ranges,npoints)
grid1_tumor, grid2_tumor = np.meshgrid(pc1_tumor, pc2_tumor)
tot_density = density_tumor - density_normal/1.5

colormap='pink_r'#'autumn_r'#'pink_r'#'gist_heat_r'#'OrRd'#'RdYlGn_r'#'RdYlBu_r'#'Oranges'#'YlOrBr'#'Greys'
npoints=25
#axs4g[0,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
#axs4g[0,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
#axs4g[1,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
#axs4g[1,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)

colormap='RdBu_r'#'seismic'#'coolwarm'
npoints=25
axs4g[0,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4g[0,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4g[1,0].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)
axs4g[1,1].contourf(grid_pc1, grid_pc2, tot_density, npoints, cmap=colormap)

colormap1='Blues'
colormap2='gist_heat_r'
al1=0.9;al2=0.55
aa=False
#axs4g[0,0].contourf(pc1_normal, pc2_normal, density_normal, npoints, cmap=colormap1, alpha=al1, antialiased=aa)
#axs4g[0,0].contourf(pc1_tumor, pc2_tumor, density_tumor, npoints, cmap=colormap2, alpha=al2, antialiased=aa)
#axs4g[0,1].contourf(pc1_normal, pc2_normal, density_normal, npoints, cmap=colormap1, alpha=al1, antialiased=aa)
#axs4g[0,1].contourf(pc1_tumor, pc2_tumor, density_tumor, npoints, cmap=colormap2, alpha=al2, antialiased=aa)
#axs4g[1,0].contourf(pc1_normal, pc2_normal, density_normal, npoints, cmap=colormap1, alpha=al1, antialiased=aa)
#axs4g[1,0].contourf(pc1_tumor, pc2_tumor, density_tumor, npoints, cmap=colormap2, alpha=al2, antialiased=aa)
#axs4g[1,1].contourf(pc1_normal, pc2_normal, density_normal, npoints, cmap=colormap1, alpha=al1, antialiased=aa)
#axs4g[1,1].contourf(pc1_tumor, pc2_tumor, density_tumor, npoints, cmap=colormap2, alpha=al2, antialiased=aa)


axs4g[0,0].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4g[0,0].scatter(sign1*pc_data[ind_stageI,0],sign2*pc_data[ind_stageI,1],label=tissue_id+', stage I', marker="s",c='r',alpha=al)
axs4g[0,1].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4g[0,1].scatter(sign1*pc_data[ind_stageII,0],sign2*pc_data[ind_stageII,1],label=tissue_id+', stage II', marker="s",c='r',alpha=al)
axs4g[1,0].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4g[1,0].scatter(sign1*pc_data[ind_stageIII,0],sign2*pc_data[ind_stageIII,1],label=tissue_id+', stage III', marker="s",c='r',alpha=al)
axs4g[1,1].scatter(sign1*pc_data[ind_normal,0],sign2*pc_data[ind_normal,1],label=tissue_id+', normal', marker="o",c='b',alpha=al)
axs4g[1,1].scatter(sign1*pc_data[ind_stageIV,0],sign2*pc_data[ind_stageIV,1],label=tissue_id+', stage IV', marker="s",c='r',alpha=al)

axs4g[1,0].set_xlabel('PC1'), axs4g[1,1].set_xlabel('PC1')
axs4g[0,0].set_ylabel('PC2'), axs4g[1,0].set_ylabel('PC2')
axs4g[0,0].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4g[0,1].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4g[1,0].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4g[1,1].legend(bbox_to_anchor=(0., 1.), edgecolor='k',loc='upper left',borderaxespad=0.2)
axs4g[0,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4g[0,0].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4g[0,1].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4g[0,1].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4g[1,0].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4g[1,0].yaxis.set_minor_locator(AutoMinorLocator(2))
axs4g[1,1].xaxis.set_minor_locator(AutoMinorLocator(2)),axs4g[1,1].yaxis.set_minor_locator(AutoMinorLocator(2))
#fig4g.colorbar()
fig4g.tight_layout()
fig4g.savefig(tissue_id + "_g_stages_fig.pdf")

print("Done!\n")
print("Check out the outputs!\n")

