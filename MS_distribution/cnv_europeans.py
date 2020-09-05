# -*- coding: utf-8 -*-
"""
At the moment, this script file reproduce some of the 
figures present in paper doi:10.1038/nature24287
refered to fitness and mutation data of the LTEE.

Simulations of the experimental data will be release
in further versions.

For more details see author(s):
DALV, AG
"""

import numpy as np
import matplotlib.pyplot as plt
from pylab import *


# Names ----------------------------------------------------------

datapath="./" 
file_cnv = "CNV.bed"


# Functions ------------------------------------------------------

def Read_Two_Column_From_File(file_name,lines,c1,c2):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        x = []
        y = []
        for line in data:
            p = line.split()
            x.append(int(p[c1-1]))
            y.append(int(p[c2-1]))

    return x, y

def Difference_Dist(x0,xf):
    if(len(x0)==len(xf)):
        diff = []
        pcum = []
        for i in range(len(x0)):
            if(xf[i]>x0[i]):
                diff.append(xf[i]-x0[i])
                pcum.append(len(x0)-i)
            else:
                print('Found not positive different at position ',i)
        diff=sort(diff)
    else:
        print("Error: length of input arrays missmatch")

    return diff, pcum


# Filling arrays ------------------------------------------------------


start_hg18, end_hg18 = Read_Two_Column_From_File(datapath+file_cnv,1,2,3)
print(start_hg18[1],end_hg18[1])
print(-start_hg18[1]+end_hg18[1])
cnv_size, cnv_cprob = Difference_Dist(start_hg18,end_hg18)


# Plots ---------------------------------------------------------------

# plotting Ws ----------------------------
plt.loglog(cnv_size, cnv_cprob, 
label = 'Cumulantive Probability',
#color='green', 
#linestyle='dashed', 
#linewidth = 3, 
marker='o' 
#, markerfacecolor='blue', markersize=12
)  

plt.xlabel('CNV size') 
plt.ylabel('Cumulant Number of CNV in European ancestry subjects')
plt.tight_layout()
plt.legend()
plt.savefig('CNV_distribution_fig.png')
#plt.show() 

