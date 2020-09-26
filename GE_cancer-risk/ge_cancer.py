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
lifetime=85
tissue_id = 'COAD'
tissues_PCA = ['BRCA','COAD','ESCA','HNSC','LIHC','LUAD','PRAD','THCA']
tissues_CLRm = ['Breast','Colorectal','Esophageal','Head_Neck','Hepatocellular','Lung','Prostate','Thyroid_Papillary']
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
        name.append(line.split()[0])
        risk.append(line.split()[-6])
        Nsc.append(line.split()[-4])
        msc.append(line.split()[-3])
    name = np.array(name)
    risk = np.array(risk, np.float)
    Nsc = np.array(Nsc, np.float)
    msc = np.array(msc, np.float)
    print(msc)
    print(Nsc)
  
    return(name,risk,Nsc,msc)

def read_CLRr(wb_name):
    name = []    
    risk = []
    error = []
    clr_wb = xl.open_workbook(wb_name)
    clr_ws = clr_wb.sheet_by_index(0)
    for i in range(2,clr_ws.ncols-5):
        name.append(clr_ws.cell_value(0,i))
        print(clr_ws.cell_value(0,i))
        risk.append(clr_ws.cell_value(-2,i))
        error.append(clr_ws.cell_value(-1,i))
    name = np.array(name)
    risk = np.array(risk, np.float)
    error = np.array(risk, np.float)
  
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
        Rn.append(np.std(pc_data[ind_normal]))
        Xt.append(abs(np.mean(pc_data[ind_tumor]) - Xn))
        Rt.append(np.std(pc_data[ind_tumor]- Xn))  
    Xt = np.array(Xt, dtype="f")      
    Rn = np.array(Rn, dtype="f")   
    Rt = np.array(Rt, dtype="f")   

    return Xt,Rn,Rt#,D


#def find_indexes()

 #   return indexes

# Reading and processing all the data ---------------------------------------

D = read_D(D_file)
print(D)
nameCLRm,risk_m,Nsc,msc = read_CLRm(CLR_path+CLRm_file)
print(nameCLRm)
nameCLRr,risk_r,error_r = read_CLRr(CLR_path+CLRr_file)
t0 = np.log2(Nsc)
t = t0 + msc*80#lifetime #total number of stem cell divisions
risk_m = risk_m/Nsc #risk per stem cell

aref = 2e-14          #reference value
ERS = risk_m/(aref*t) #extra risk score

#working with PCA data
pc_data,ind_normal,ind_tumor = read_PC(sample_path,PC_path,tissue_id,Npc)
pc,mfitness = fitness_dist(-pc_data[:,0],ind_normal,ind_tumor,1.,1.5,5,16)
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_PCA)
R=Xt-(Rt+Rn)

indexes_m= []
for t_clr in tissues_CLRm:
    for i in range(len(nameCLRm)):
        if t_clr==nameCLRm[i]: indexes_m.append(i)
print(nameCLRm[indexes_m])
indexes_r= []
for t_clr in tissues_CLRr:
    for i in range(len(nameCLRr)):
        if t_clr in nameCLRr[i]: indexes_r.append(i)
print(nameCLRr[indexes_r])

log_risk = np.log(risk_r[indexes_r]/Nsc[indexes_m])
print(log_risk)
print(np.log(D*t[indexes_r]/R))


# Exporting readable data ------------------------------------------------------

sl = max(np.char.str_len(nameCLRm)) #maximum number of characters of nameCLRm

with open('ERS_dat.dat','w+') as savefile:
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

with open('GE_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of GE for each tissue.\n\n')
    savefile.write('Tissue\tXt\t\tRn\t\tRt\t\tR\t\tD\t\tNsc\t\tmsc\t\trisk\n')
    for i in range(len(tissues_PCA)):
        savefile.write(' '+tissues_PCA[i]+'\t')
        savefile.write('%f\t' % Xt[i])
        savefile.write('%f\t' % Rn[i])
        savefile.write('%f\t' % Rt[i])
        savefile.write('%f\t' % R[i])
        savefile.write('%f\t' % D[i])
        savefile.write('%f\t' % Nsc[indexes_m[i]])
        savefile.write('%f\t' % msc[indexes_m[i]])
        savefile.write('%f\n' % risk_r[indexes_r[i]])
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

plt.loglog(t,risk_m,linestyle='none',marker='o')
for i in range(len(t)):
    plt.annotate(nameCLRm[i], (t[i], risk_m[i]))

plt.xlabel('t0+msc*80yrs')
plt.ylabel('Lifetime risk/Nsc')
plt.tight_layout()
plt.savefig('lifetime-risk_fig.pdf')
plt.clf() # clear the plot

   
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
plt.clf() # clear the plot

plt.plot(pc,mfitness,label='fitness',c='black')
plt.xlabel('PC1')
plt.ylabel('-fitness')
plt.title(tissue_id)
plt.legend()
plt.tight_layout()
plt.savefig('ge'+tissue_id+'_fitness_fig.pdf')
plt.clf() # clear the plot

# plotting Brownian and Levy Correlation ----------------------------

#Calculate x for Brownian oscilations plot
x_b = np.array(-2*(D*t[indexes_m]**0.5/R)**-2 + np.log(D*t[indexes_m]**0.5/R), dtype="f")  

#Calculate x for large Levy jumps
x_l = np.array(np.log(D*t[indexes_m]/R), dtype="f")  

#Prepare scatter plot for Brownian oscilations
plt.scatter(x_b, log_risk, label = 'Data')

#Linear fit for Brownian oscilations
m, b = np.polyfit(x_b, log_risk, 1)

p2 = plt.plot([np.min(x_b),np.max(x_b)], [m*np.min(x_b),m*np.max(x_b)]+b,
color='r', 
label = 'Linear fit',
linestyle='dashed')

#Set plot parameters
plt.xlabel('-2(Dt^0.5/R)^-2 + ln(Dt**0.5/R)')
plt.ylabel('ln(Risk/Nsc)')
plt.tight_layout()
plt.legend(loc=3)
plt.savefig('brownian_fig.pdf')
#plt.show()
plt.clf() # clear the plot

#Prepare scatter plot for large Levy jumps
plt.scatter(x_l, log_risk, label = 'Data')
for i in range(len(x_l)):
    plt.annotate(tissues_PCA[i], (x_l[i], log_risk[i]))

#Linear fit for large Levy jumps
m = 1
n = np.mean(log_risk - m*x_l)
p2 =plt.plot([np.min(x_l),np.max(x_l)], m*[np.min(x_l),np.max(x_l)]+n,
color='r',
label ='-22.91+x',
linestyle='dashed')

#Set plot parameters
plt.xlabel('-2Dt/Rn^-2 + ln(Dt/Rn)')
plt.ylabel('ln(Risk/Nsc)')
plt.tight_layout()
plt.legend(loc=4)
plt.savefig('levy_fig.pdf')

