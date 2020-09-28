"""
This script file processes data on cancer risk and the
results of the PCA technique applied to TCGA data of a 
set of tissues. The data surces are mentioned in the 
imported files, while the procedures and main results 
are described in paper arxiv:1507.08232v4

For more details on this file see author(s):
RHP, DALV
"""

import numpy as np
import xlrd as xl
import matplotlib as mpl
import matplotlib.pyplot as plt


# General variables -----------------------------------------------

Npc = 2
lifetime=80
tissue_id = 'COAD'
tissues_PCA = ['BRCA','COAD','ESCA','HNSC','LIHC','LUAD','PRAD','THCA']
tissues_CLRm = ['Breast','Colorectal','Esophageal','Head&Neck','Hepatocellular','Lung','Prostate','Thyroid_Papillary']
tissues_CLRr = ['breast','colorectal','esophageal','head & neck','hepatocellular','lung','prostate','thyroid follicular']

sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'
CLR_path = '../databases_external/CLR/'
CLRm_file = '1260825_Table1.dat'
CLRr_file = 'aaf9011_Table_y'+str(lifetime)+'.xlsx'
D_file = 'tissues_D-data.txt'


# Functions ------------------------------------------------------

def read_D(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()

    D = []    
    for line in lines[1:]:
        D.append(line.split()[-1])
    D = np.array(D, np.float)
  
    return(D)


def read_CLRm(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()

    name = []    
    risk = []
    Nsc = []
    msc = []
    for line in lines[16:]:
        print(line.split()[0])
        name.append(line.split()[0])
        risk.append(line.split()[-6])
        Nsc.append(line.split()[-4])
        msc.append(line.split()[-3])
    name = np.array(name)
    risk = np.array(risk, np.float)
    Nsc = np.array(Nsc, np.float)
    msc = np.array(msc, np.float)
  
    return(name,risk,Nsc,msc)

def read_CLRr(wb_name):
    name = []    
    risk = []
    error = []
    clr_wb = xl.open_workbook(wb_name)
    clr_ws = clr_wb.sheet_by_index(0)
    for i in range(2,clr_ws.ncols-5):
        print(clr_ws.cell_value(0,i))
        name.append(clr_ws.cell_value(0,i))
        risk.append(clr_ws.cell_value(-2,i))
        error.append(clr_ws.cell_value(-1,i))
    name = np.array(name)
    risk = np.array(risk, np.float)
    error = np.array(error, np.float)
  
    return(name,risk,error)
          
def read_PC(sample_path,PC_path,tissue_id,Npc):
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
        for j in range(Npc):    
            pc_data[i].append(pc_ws.cell_value(i,j))

    pc_data = np.array(pc_data, dtype="f")

    return pc_data, ind_normal, ind_tumor


def fitness_dist(pc1,ind_n,ind_t,fitN,fitT,binsN,binsT):

    histN,bin_edgesN = np.histogram(pc1[ind_n],bins=binsN)
    histT,bin_edgesT = np.histogram(pc1[ind_t],bins=binsT)
    pc = np.concatenate((bin_edgesN[0:-1],bin_edgesT[0:-1]), axis=None)
    mfitness = np.concatenate((-fitN*histN/max(histN),-fitT*histT/max(histT)))

    return pc,mfitness


def compute_GEdistances(sample_path,PC_path,tissues):
    Xt = []
    Rn = []
    Rt = []
    #D = []
    for t_id in tissues:
        pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,t_id,1)
        Xn = np.mean(pc_data[ind_normal])
        Rn.append(np.std(pc_data[ind_normal],ddof=1))
        Xt.append(abs(np.mean(pc_data[ind_tumor]) - Xn))
        Rt.append(np.std(pc_data[ind_tumor]-Xn,ddof=1))  
    Xt = np.array(Xt, dtype="f")      
    Rn = np.array(Rn, dtype="f")   
    Rt = np.array(Rt, dtype="f")   

    return Xt,Rn,Rt#,D


# Reading and processing all the data ---------------------------------------

D = read_D(D_file)
print("Loading lifetime variables for tissues:")
nameCLRm,risk_m,Nsc,msc = read_CLRm(CLR_path+CLRm_file)
print("Done!\nLoading cancer risks for tissues:")
nameCLRr,risk_r,error_r = read_CLRr(CLR_path+CLRr_file)
t0 = np.log2(Nsc)
t_r = t0 + msc*lifetime #total number of stem cell divisions
t_m = t0 + msc*80
risk_m = risk_m/Nsc     #risk per stem cell
aref = 2e-14            #reference value
ERS = risk_m/(aref*t_m) #extra risk score

print("Done!\nLoading Principal Componets for:") #working with PCA data
pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissue_id,Npc)
pc,mfitness = fitness_dist(-pc_data[:,0],ind_normal,ind_tumor,1.,1.5,5,16)
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_PCA)
R=Xt-(Rt+Rn)

indexes_m= []
for t_clr in tissues_CLRm:
    for i in range(len(nameCLRm)):
        if t_clr==nameCLRm[i]: indexes_m.append(i)
indexes_r= []
for t_clr in tissues_CLRr:
    for i in range(len(nameCLRr)):
        if t_clr in nameCLRr[i]: indexes_r.append(i)

log_risk = np.log(risk_r[indexes_r]/Nsc[indexes_m])

x_brownian=np.array(-2*(D*t_m[indexes_m]**0.5/R)**-2+np.log(D*t_m[indexes_m]**0.5/R), dtype="f")  
x_levy = np.array(np.log(D*t_r[indexes_m]/R), dtype="f")  
m_brownian, b_brownian=np.polyfit(x_brownian,log_risk,1) #Linear fit for Brownian
m_levy=1; n_levy=np.mean(log_risk-m_levy*x_levy)         #Linear fit for Levy jumps


# Exporting readable data ------------------------------------------------------

sl = max(np.char.str_len(nameCLRm)) #maximum number of characters of nameCLRm

with open('geERS_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of ERS for each tissue.\n\n')
    savefile.write('#Tissue')
    for i in range(sl-6):
        savefile.write(' ')
    savefile.write('\tERS\n')
    for i in range(len(nameCLRm)):
        savefile.write(" %s" %nameCLRm[i])
        for k in range(sl-len(nameCLRm[i])):
            savefile.write(' ')
        savefile.write("\t%f\n" %ERS[i])
    savefile.write('\n')
    savefile.close()

with open('geData_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of GE for each tissue.\n\n')
    savefile.write('Tissue\tXt\t\tRn\t\tRt\t\tR\t\tD\t\t<risk>\t\t<dev>\n')
    for i in range(len(tissues_PCA)):
        savefile.write(' '+tissues_PCA[i]+'\t')
        savefile.write('%f\t' % Xt[i])
        savefile.write('%f\t' % Rn[i])
        savefile.write('%f\t' % Rt[i])
        savefile.write('%f\t' % R[i])
        savefile.write('%f\t' % D[i])
        savefile.write('%f\t' % risk_r[indexes_r[i]])
        savefile.write('%f\n' % error_r[indexes_r[i]])
    savefile.write('\n')
    savefile.close()


# Plots ---------------------------------------------------------------

# plotting time vs. risk ----------------------------
r11 = 1e-12
r21 = 2e-13
r12 = 1e-9
r22 = 2e-10
x1 = 1e1
x2 = 1e4

plt.loglog([x1, x2],[r11, r12],color='r',linestyle = 'dashed')
plt.loglog([x1, x2],[r21, r22],color='r',linestyle = 'dashed')

plt.loglog(t_m,risk_m,linestyle='none',marker='o')
for i in range(len(t_m)):
    plt.annotate(nameCLRm[i],(t_m[i],risk_m[i]), ha='center')

plt.xlabel('t0+msc*'+str(80)+'yrs')
plt.ylabel('Lifetime risk/Nsc')
plt.tight_layout()
plt.savefig('risk-lifetime_fig.pdf')
plt.clf()
   
# plotting PC1 vs. PC2 ----------------------------
plt.scatter(-pc_data[ind_normal,0],pc_data[ind_normal,1],label='normal',c='b',s=15)
plt.scatter(-pc_data[ind_tumor,0],pc_data[ind_tumor,1],label='tumor',c='r',s=15)
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.title(tissue_id)
plt.legend()
plt.tight_layout()
plt.grid(True)
plt.savefig('ge'+tissue_id+'_pc1-2_fig.pdf')
plt.clf()

plt.plot(pc,mfitness,label='fitness',c='black')
plt.xlabel('PC1')
plt.ylabel('-fitness')
plt.title(tissue_id)
plt.legend()
plt.tight_layout()
plt.savefig('ge'+tissue_id+'_fitness_fig.pdf')
plt.clf()

# plotting Brownian and Levy Correlation ----------------------------
plt.scatter(x_brownian, log_risk, label = 'Data')
for i in range(len(x_brownian)):
    plt.annotate(tissues_PCA[i], (x_brownian[i], log_risk[i]), ha='center')

plt.plot([np.min(x_brownian),np.max(x_brownian)], 
[m_brownian*np.min(x_brownian),m_brownian*np.max(x_brownian)]+b_brownian,
color='r', label='Linear fit', linestyle='dashed')

plt.xlabel('-2(Dt^0.5/R)^-2 + ln(Dt^0.5/R)')
plt.ylabel('ln(Risk/Nsc)')
plt.tight_layout()
plt.legend(loc=3)
plt.savefig('risk-brownian_fig.pdf')
plt.clf()

plt.scatter(x_levy, log_risk, label = 'Data')
for i in range(len(x_levy)):
    plt.annotate(tissues_PCA[i], (x_levy[i], log_risk[i]), ha='center')
plt.plot([np.min(x_levy),np.max(x_levy)], m_levy*[np.min(x_levy),np.max(x_levy)]+n_levy,
color='r', label=str(n_levy)+'+'+str(m_levy)+'x', linestyle='dashed')

plt.xlabel('ln(Dt/R)')
plt.ylabel('ln(Risk/Nsc)')
plt.tight_layout()
plt.legend(loc=4)
plt.savefig('risk-levy_fig.pdf')

print("Done!")
