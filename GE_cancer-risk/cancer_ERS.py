"""
This script file processes data on cancer risk for a set
of tissues. A reference lower limit value is estimated
and the Extra Risk Score is calculated for all tissues.
The procedure and main results are described in paper 
arxiv:1507.08232v4

For more details on this file see author(s):
RHP, DALV
"""

import numpy as np
import matplotlib.pyplot as plt

# Names ----------------------------------------------------------

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
    risk = np.array(risk, np.float)
    Nsc = np.array(Nsc, np.float)
    msc = np.array(msc, np.float)

  
    return(name,risk,Nsc,msc)


# Reading and processing the data ------------------------------------------

name,risk,Nsc,msc = read_CLRdata(CLR_path+CLR_file)
t0 = np.log2(Nsc)
t = t0 + msc*80 #total number of stem cell divisions
risk = risk/Nsc #risk per stem cell

aref = 2e-14             #reference value
ERS = risk/(aref*(t0+t)) #extra risk score
    

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
plt.savefig('lifetimerisk_fig.pdf')
   

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

