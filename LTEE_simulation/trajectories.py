"""
At the moment, this script file reproduce some of the 
figures present in paper doi:10.1038/nature24287
refered to fitness and mutation data of the LTEE.

Simulations of the experimental data will be release
in further versions.

For more details on this file see author(s):
DALV
"""

import numpy as np
import matplotlib.pyplot as plt


# Names ----------------------------------------------------------

datapath="../databases_external/LTEE/trajectories/" 
file_Ws = "Ws.dat"
file_Ms = "Ms.dat"


# Functions ------------------------------------------------------

def Read_Trajectories_File(file_name,lines):
    with open(file_name, 'r') as rfile:
        data = rfile.readlines()[lines:]
        popul = []
        num = []
        gen = []
        traj = []
        ipop = 0 
        gen.append([])
        traj.append([])
        for line in data:
            if line.startswith("p"):
            	p = line.split()
            	num.append(0)
            	popul.append(p[1])
            elif line =='\n': 
            	ipop = ipop + 1
            	gen.append([])
            	traj.append([])
            else:
            	p = line.split()
            	num[ipop] = num[ipop] + 1
            	gen[ipop].append(int(p[0]))
            	traj[ipop].append(float(p[1]))

        #for i in range(len(size)):
        #    size[i] = sort(size[i])

    return popul, num, gen, traj


# Filling arrays ------------------------------------------------------


populWs, numWs, genWs, trajWs = Read_Trajectories_File(datapath+file_Ws,2)
populMs, numMs, genMs, trajMs = Read_Trajectories_File(datapath+file_Ms,2)


# Plots ---------------------------------------------------------------

# plotting Ws ----------------------------
for i in range(len(populWs)):
    plt.plot(genWs[i], trajWs[i], 
    label = populWs[i],
    #color='green', 
    #linestyle='dashed', 
    #linewidth = 3, 
    marker='o' 
    #, markerfacecolor='blue', markersize=12
    )  

plt.xlabel('generations') 
plt.ylabel('Fitness')
plt.tight_layout()
plt.legend()
plt.savefig('Ws_fig.png')
plt.clf()


# plotting Ms ----------------------------
 
for i in range(len(populMs)):
    plt.plot(genMs[i], trajMs[i], 
    label = populMs[i],
    #color='green', 
    #linestyle='dashed', 
    #linewidth = 3, 
    marker='o' 
    #, markerfacecolor='blue', markersize=12
    )  

plt.xlabel('generations') 
plt.ylabel('Mutations')
plt.tight_layout()
plt.legend()
plt.savefig('Ms_fig.png')


# Quick Note ----------------------------
print("\nAt the moment we are plotting just the \nexperimental data we want to reproduce")
print("Simulations coming soon!\n")

