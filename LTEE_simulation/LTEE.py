# -*- coding: utf-8 -*-
"""
This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

from pylab import *
#from netCDF4 import Dataset

datapath="../external_databases/LTEE/" 
file_mutF = "Ara-1_mut-fixed.dat"
file_SNP = "Ara-1_mut-SNP.dat"

def Read_Two_Column_File(file_name,lines):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        x = []
        y = []
        for line in data:
            p = line.split()
            x.append(int(p[0]))
            y.append(int(p[1]))

    return x, y

def Read_SNP_File(file_name,lines):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        gen = []
        fr = []
        igen = 0 
        fr.append([])
        for line in data:
            if line.startswith("f"):
            	p = line.split()
            	gen.append(int(p[3]))
            elif line =='\n': 
            	igen = igen + 1
            	fr.append([])
            else:
            	fr[igen].append(float(line))

    return gen, fr

def Compute_MutT(genF, mutF, genSNP, fr):
    if genSNP != genF:
    	print("The sampling of the SNP doesn't match the sampling of the fixed mutations")
    mutT = []
    for i in range(len(genF)):
    	mutT.append( mutF[i] + sum(fr[i]) ) # similar to paper: Mutations as levy flights

    return mutT

#------------------------------------------------------


genF, mutF = Read_Two_Column_File(datapath+file_mutF,4)
genSNP, fr = Read_SNP_File(datapath+file_SNP,3)
mutT = Compute_MutT(genF, mutF, genSNP, fr)



# Plots ..............................................
# plotting the points 
plt.plot(genF[0:5], mutF[0:5], 
label = "fixed",
#color='green', 
linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.plot(genF[0:5], mutT[0:5], 
label = "total",
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)

plt.xlabel('gen') 
plt.ylabel('Number of SPMs')
plt.tight_layout()
plt.legend()
plt.xlim(0, 20250)
plt.savefig('mut_20K_fig.pdf')
#plt.show() 

# clear the plot
plt.clf()

# plotting up to 40K 
plt.axvline(x=27000,color='gray',linestyle='dotted')

plt.plot(genF, mutF, 
label = "fixed",
#color='green', 
linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.plot(genF, mutT, 
label = "total",
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)

plt.xlabel('gen') 
plt.ylabel('Number of SPMs')
plt.tight_layout()
plt.legend()
plt.savefig('mut_40K_fig.pdf')
#plt.show() 
