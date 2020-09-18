# -*- coding: utf-8 -*-
"""
This script file processes available CNV data 
of over 100,000 European ancestry subjects.
The output figure is similar to the one in paper
arXiv:1605.09697

For more details on this file see author(s):
DALV
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
from os.path import *
from pylab import *


# Names ----------------------------------------------------------

datapath="../databases_external/CNV/" 
file_cnv = "eur_CNV.bed"


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
                pcum.append((len(x0)-i)/len(x0))
            else:
                print('Found not positive different at position ',i)
        diff=sort(diff)
    else:
        print("Error: length of input arrays missmatch")

    return diff, pcum


# Manipulating the data -----------------------------------------------

print("Processing CNV data of over 100,000 European ancestry subjects ...")
if(not isfile(datapath+file_cnv)):
    print("CNV databases not found, please read", "'"+ datapath + "info.txt'")
    sys.exit()

start_hg18, end_hg18 = Read_Two_Column_From_File(datapath+file_cnv,1,2,3)
print(len(start_hg18),"values were loaded!")
print("10 secs more to compute and produce the plot")
cnv_size, cnv_cprob = Difference_Dist(start_hg18,end_hg18)

# Plots ---------------------------------------------------------------

# plotting dist ----------------------------
plt.loglog(cnv_size, cnv_cprob, 
label = 'Cumulantive Number',
color='#084594',
#mfc='none',
linestyle='none', 
#linewidth = 3, 
marker='o' 
, markerfacecolor='none', markersize=2
)  

plt.xlabel('CNV size') 
plt.ylabel('Cumulative Number of CNV in European ancestry subjects')
plt.tight_layout()
plt.legend(loc='upper right')
plt.savefig('eur_CNV_distribution_fig.png')
#plt.show()

print("Job Done!") 

