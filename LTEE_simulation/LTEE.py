# -*- coding: utf-8 -*-
"""
This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt

from pylab import *
#from netCDF4 import Dataset

datapath="../external_databases/LTEE/" 
file_A1_mutF = "Ara-1_mut-SP_fixed.dat"
file_A1_SNP = "Ara-1_mut-SNP.dat"
file_A1_LR = "Ara-1_mut-LR.dat"

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

def Read_A1_LR_File(file_name,lines):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        gen = []
        num = []
        size = []
        igen = 0 
        size.append([])
        for line in data:
            if line.startswith("g"):
            	p = line.split()
            	num.append(0)
            	gen.append(int(p[1]))
            elif line =='\n': 
            	igen = igen + 1
            	size.append([])
            else:
            	num[igen] = num[igen] + 1
            	size[igen].append(int(line))

    return gen, num, size

def Compute_MutT(genF, mutF, genSNP, fr):
    if genSNP != genF:
    	print("The sampling of the SNP doesn't match the sampling of the fixed mutations")
    mutT = []
    for i in range(len(genF)):
    	mutT.append( mutF[i] + sum(fr[i]) ) # similar to paper: Mutations as levy flights

    return mutT

#------------------------------------------------------


genF, mutF = Read_Two_Column_File(datapath+file_A1_mutF,4)
genSNP, fr = Read_SNP_File(datapath+file_A1_SNP,3)
mutT = Compute_MutT(genF, mutF, genSNP, fr)
genLR, numLR, sizeLR = Read_A1_LR_File(datapath+file_A1_LR,3)


#------------------------------------------------------

with open("Ara-1_mut-SP_total_dat.dat", "w+") as file:
    file.write("Number of total mutations in mixed-population samples\n\n")
    file.write(" gen\tnum\n")
    for x in zip(genF, mutT):
        file.write(" {0}\t{1}\n".format(*x))


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
plt.savefig('mut-SP_20K_fig.pdf')
#plt.show() 

# clear the plot
plt.clf()

# plotting SP up to 40K 
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
plt.savefig('mut-SP_40K_fig.pdf')
#plt.show() 


# clear the plot
plt.clf()

# plotting LR up to 50K 

plt.plot(genLR, numLR, 
label = "Ara-1",
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.xlabel('gen') 
plt.ylabel('Number of Rearrangements')
plt.tight_layout()
plt.legend()
plt.savefig('mut-LR_50K_fig.pdf')
#plt.show() 

