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
file_Aall_LR = "Ara-all_mut-LR.dat"

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

def Read_Aall_LR_File(file_name,lines):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        popul = []
        num = []
        size = []
        ipop = 0 
        size.append([])
        for line in data:
            if line.startswith("p"):
            	p = line.split()
            	num.append(0)
            	popul.append(p[1])
            elif line =='\n': 
            	ipop = ipop + 1
            	size.append([])
            else:
            	num[ipop] = num[ipop] + 1
            	size[ipop].append(int(line))

        for i in range(len(size)):
            size[i] = sort(size[i])

    return popul, num, size


def Compute_MutT(genF, mutF, genSNP, fr):
    if genSNP != genF:
    	print("Datadases: The sampling of the SNP doesn't match the sampling of the fixed mutations")
    mutT = []
    for i in range(len(genF)):
    	mutT.append( mutF[i] + sum(fr[i]) ) # similar to paper: Mutations as levy flights

    return mutT

def Fill_SizeLR(sizeLR):

    def Included_In(value, array):
        included = False
        for i in range(len(array)):
            if value == array[i]:
            	included = True
            	break
        return included

    sizeLRu = []
    cnumLRu = []
    sizeLRu.append(sizeLR[0][0])
    for i in range(len(sizeLR)):
        for j in range(len(sizeLR[i])):
            if not Included_In(sizeLR[i][j], sizeLRu):
            	sizeLRu.append(sizeLR[i][j])
    sizeLRu = sort(sizeLRu)
    for i in range(len(sizeLRu)):
    	cnumLRu.append( len(sizeLRu) - i )

    return sizeLRu, cnumLRu

# Filling arrays ------------------------------------------------------


genF, mutF = Read_Two_Column_File(datapath+file_A1_mutF,4)
genSNP, fr = Read_SNP_File(datapath+file_A1_SNP,3)
mutT = Compute_MutT(genF, mutF, genSNP, fr)
genLR, numLR, sizeLR = Read_A1_LR_File(datapath+file_A1_LR,3)
sizeLRu, cnumLRu = Fill_SizeLR(sizeLR)
populN, numLRall, sizeLRall = Read_Aall_LR_File(datapath+file_Aall_LR,3)
sizeLRjoined, cnumLRjoined = Fill_SizeLR(sizeLRall)
for i in range(len(cnumLRjoined)):
    cnumLRjoined[i] = cnumLRjoined[i]/len(populN)

#print(sizeLRjoined)
#print(cnumLRjoined)

#------------------------------------------------------

with open("Ara-1_mut-SP_total_dat.dat", "w+") as file:
    file.write("Number of total mutations in mixed-population samples\n\n")
    file.write(" gen\tnum\n")
    for x in zip(genF, mutT):
        file.write(" {0}\t{1}\n".format(*x))


#.............................................#
#                                             #
#                  Plots                      #
#                                             #
#.............................................#

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

# clear the plot
plt.clf()

# plotting LR sizes per generation

for i in range(len(numLR)):
    sizeLRi = sizeLR[i]
    pnum = []
    for j in range(numLR[i]):
        pnum.append((numLR[i]-j)/genLR[i])

    plt.loglog(sizeLRi, pnum, 
    label = "gen " + str(genLR[i]),
    #color='green', 
    #linestyle='dashed', 
    #linewidth = 3, 
    marker='o' 
    #, markerfacecolor='blue', markersize=12
    )  

plt.xlabel('size') 
plt.ylabel('Ara-1: Cumulant Number of Rearrangements / Generation')
plt.tight_layout()
plt.legend()
plt.savefig('mut-LR_size_fig.pdf')
#plt.show() 

# clear the plot
plt.clf()

# plotting LR sizes
plt.loglog(sizeLRu, cnumLRu, 
label = "Ara-1 in 2-50K gen",
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.xlabel('size') 
plt.ylabel('Ara-1: Cumulant Number of All Rearrangements')
plt.tight_layout()
plt.legend()
plt.savefig('mut-LR_sizeU_fig.pdf')
#plt.show() 


# clear the plot
plt.clf()

# plotting LR Ara-all each sizes at 40000 generations
for i in range(len(numLRall)):
    pnum = []
    for j in range(numLRall[i]):
        pnum.append(numLRall[i]-j)

    plt.loglog(sizeLRall[i], pnum, 
    label = populN[i],
    #color='green', 
    #linestyle='dashed', 
    #linewidth = 3, 
    marker='o' 
    #, markerfacecolor='blue', markersize=12
    )  

plt.xlabel('size') 
plt.ylabel('Ara: Cumulant Number of Rearrangements at 40K Generations')
plt.tight_layout()
plt.legend()
plt.savefig('mut-LRall_size_fig.pdf')
#plt.show() 

# clear the plot
plt.clf()

# plotting LR Ara-all joined sizes at 40000 generations
plt.loglog(sizeLRjoined, cnumLRjoined, 
label = "Ara-all at 40K gen",
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.xlabel('size') 
plt.ylabel('Ara-all: Cumulant Number of All Rearrangements')
plt.tight_layout()
plt.legend()
plt.savefig('mut-LRallj_size_fig.pdf')
#plt.show() 


# plotting LR Ara-all each sizes at 40000 generations
for i in range(len(numLRall)):
    pnum = []
    for j in range(numLRall[i]):
        pnum.append(numLRall[i]-j)

    plt.loglog(sizeLRall[i], pnum, 
    label = populN[i],
    #color='green', 
    #linestyle='dashed', 
    #linewidth = 3, 
    marker='o' 
    #, markerfacecolor='blue', markersize=12
    )  

plt.xlabel('size') 
plt.ylabel('Ara: Cumulant Number of Rearrangements at 40K Generations')
plt.tight_layout()
plt.legend()
#plt.savefig('mut-LRall_size_fig.pdf')
plt.show() 

