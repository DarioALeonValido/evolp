import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('default')

#Data import

def read_data(file_name = 'cancer_risk.dat'):
    file = open(file_name, 'r')
    lines = file.readlines()
    
    data = []
    for line in lines[13:]:
        data.append(line.split()[-6:])
    
    data = np.array(data, dtype=np.float)

    names = []
    for line in lines[13:]:
        names.append(line.split()[0])
    
    return(data,names)

if __name__ == "__main__":
    data, names = read_data()
    risk = data[:,0]
    Nsc = data[:,2]
    msc = data[:,3]
    t0 = np.log2(Nsc)
    t = msc*80
    
#Calculate y for plot
    y = np.log10(risk/Nsc)

#Calculate x for plot
    x = np.log10(t0 + t)    


#Plot
m = 1
b1 = -13
b2 = -13.7
x1 = (1.25,4)
x1 = np.array(x1)
plt.scatter(x,y)
plt.plot(x1, x1*m + b1,color='r',linestyle = 'dashed')
plt.plot(x1, x1*m + b2,color='r',linestyle = 'dashed')
plt.xlabel('log(t0+msc80yrs)')
plt.ylabel('log(Lifetime risk/Nsc)')
plt.tight_layout()
#plt.show()
plt.savefig('lifetimerisk_fig.pdf')
   
#Calculate ERS
aref = 2e-14
ERS = risk/(Nsc*aref*(t0+t))
print('ERS: \n', ERS)    

#Save data
D = len(names)
Iter = range(D)
savefile = open('ERS_dat.dat','w')
savefile.write('In this file are listed the values of ERS for each tissue \n \n')
savefile.write('Tissue \t ERS \n')
for i in Iter:
    print(names[i], ERS[i], file = savefile)
savefile.close()
