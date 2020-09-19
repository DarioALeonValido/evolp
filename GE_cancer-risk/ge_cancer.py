"""
This script file processes data on cancer risk and the
results of the PCA technique applied to TCGA data of a 
set of tissues. The procedures and main results are 
described in paper arxiv:1507.08232v4

For more details on this file see author(s):
RHP, DALV
"""

import numpy as np
import xlrd as xl
import matplotlib.pyplot as plt

# General variables -----------------------------------------------

Npc = 2
tissue_id = 'COAD'
tissues_PCA = ['BRCA','COAD','HNSC','LIHC','LUAD','PRAD','THCA','ESCA']
tissues_CLR = ['Breast','Colorectal','Head_Neck','Hepatocellular','Lung','Prostate','Thyroid_Papillary','Esophageal']

sample_path = '../databases_external/TCGA/'
PC_path = '../databases_generated/TCGA_pca/'
CLR_path = '../databases_external/CLR/'
CLR_file = 'cancer-risk.dat'


# Functions ------------------------------------------------------

def read_CLRdata(file_name):
    file = open(file_name, 'r')
    lines = file.readlines()

    name = []    
    risk = []
    Nsc = []
    msc = []
    for line in lines[13:]:
        name.append(line.split()[0])
        risk.append(line.split()[-6])
        Nsc.append(line.split()[-4])
        msc.append(line.split()[-3])
    name = np.array(name)
    risk = np.array(risk, np.float)
    Nsc = np.array(Nsc, np.float)
    msc = np.array(msc, np.float)
  
    return(name,risk,Nsc,msc)

          
def read_PCdata(sample_path,PC_path,tissue_id,Npc):
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


def compute_GEdistances(sample_path,PC_path,tissues):
    Xt = []
    Rn = []
    Rt = []
    #D = []
    for t_id in tissues:
        pc_data,ind_normal,ind_tumor = read_PCdata(sample_path,PC_path,t_id,1)
        Xn = np.mean(pc_data[ind_normal])
        Rn.append(np.std(pc_data[ind_normal]))
        Xt.append(np.mean(pc_data[ind_tumor]) - Xn)
        Rt.append(np.std(pc_data[ind_tumor]- Xn))
  
    Xt = np.array(Xt, dtype="f")      
    Rn = np.array(Rn, dtype="f")   
    Rt = np.array(Rt, dtype="f")   

    return Xt,Rn,Rt#,D


# Reading and processing all the data ---------------------------------------

name,risk,Nsc,msc = read_CLRdata(CLR_path+CLR_file)
t0 = np.log2(Nsc)
t = t0 + msc*80 #total number of stem cell divisions
risk = risk/Nsc #risk per stem cell

aref = 2e-14             #reference value
ERS = risk/(aref*(t0+t)) #extra risk score

#working with PCA data
pc_data,ind_normal,ind_tumor = read_PCdata(sample_path,PC_path,tissue_id,Npc)
Xt,Rn,Rt = compute_GEdistances(sample_path,PC_path,tissues_PCA)


# Exporting readable data ------------------------------------------------------

sl = max(np.char.str_len(name)) #maximum number of characters of name

with open('ERS_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of ERS for each tissue.\n\n')
    savefile.write('#Tissue')
    for i in range(sl-6):
        savefile.write(' ')
    savefile.write('\tERS\n')
    for i in range(len(name)):
        savefile.write(" %s" %name[i])
        for k in range(sl-len(name[i])):
            savefile.write(' ')
        savefile.write("\t%f\n" %ERS[i])
    savefile.write('\n')
    savefile.close()

with open('GE_dat.dat','w+') as savefile:
    savefile.write('In this file are listed the values of GE for each tissue.\n\n')
    savefile.write('Tissue\tXt\t\tRn\t\tRt\n')
    for i in range(len(tissues_PCA)):
        savefile.write(' '+tissues_PCA[i]+'\t')
        savefile.write('%f\t' % Xt[i])
        savefile.write('%f\t' % Rn[i])
        savefile.write('%f\n' % Rt[i])
    savefile.write('\n')
    savefile.close()


# Plots ---------------------------------------------------------------

r11 = 1e-12
r21 = 2e-13
r12 = 1e-9
r22 = 2e-10
x1 = 1e1
x2 = 1e4

plt.loglog([x1, x2],[r11, r12],color='r',linestyle = 'dashed')
plt.loglog([x1, x2],[r21, r22],color='r',linestyle = 'dashed')

plt.loglog(t,risk,linestyle='none',marker='o')
for i in range(len(t)):
    plt.annotate(name[i], (t[i], risk[i]))

plt.xlabel('t0+msc*80yrs')
plt.ylabel('Lifetime risk/Nsc')
plt.tight_layout()
plt.savefig('lifetime-risk_fig.pdf')
plt.clf() # clear the plot
   

plt.scatter(-pc_data[ind_normal,0],pc_data[ind_normal,1],label='normal')
plt.scatter(-pc_data[ind_tumor,0],pc_data[ind_tumor,1],label='tumor')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.tight_layout()
plt.savefig('ge'+tissue_id+'_fig.pdf')
plt.clf() # clear the plot


#Import data for ploteable tissues
D = np.loadtxt('tissues_data.txt')
D = np.array(D, dtype="f")   

indexes= []
for t_clr in tissues_CLR:
    ind = np.where(name == t_clr)
    indexes.append(ind[0][0])

log_risk = np.log(risk[indexes])

#Calculate x for Brownian oscilations plot
x_b = np.array(-2*(D*t[indexes]/Rn)**(-2) + np.log(D*t[indexes]/Rn), dtype="f")  

#Calculate x for large Levy jumps
x_l = np.array(np.log(D*t[indexes]/Rn), dtype="f")  

#Prepare scatter plot for Brownian oscilations
plt.scatter(x_b, log_risk, label = 'Data')

#Linear fit for Brownian oscilations
m, b = np.polyfit(x_b, log_risk, 1)

p2 = plt.plot([np.min(x_b),np.max(x_b)], [m*np.min(x_b),m*np.max(x_b)]+b,
color='r', 
label = 'Linear fit',
linestyle='dashed')
print('Slope for brownian oscilations m =',+m)

#Set plot parameters
plt.xlabel('-2Dt/Rn^-2 + ln(Dt/Rn)')
plt.ylabel('ln(Risk/Nsc)')
plt.tight_layout()
plt.legend(loc=3)
plt.savefig('brownian_fig.pdf')
#plt.show()
plt.clf() # clear the plot

#Prepare scatter plot for large Levy jumps
#p1 = plt.scatter(x_l, y, label = 'Data')
plt.scatter(x_l, log_risk, label = 'Data')

#Linear fit for large Levy jumps
m = 1
n = np.mean(log_risk - m*x_l)
print('n=',n)
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

