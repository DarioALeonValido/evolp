import numpy as np
import matplotlib.pyplot as plt
from matplotlib import style
style.use('bmh')
          
#Impot data for ploteable tissues
x1, Rn, Rt, D, Nsc, msc, risk = np.loadtxt('tissues_data.txt', delimiter=',', unpack=True)

#Calculate y for plot
y = np.log(risk/Nsc)

#Calculate t0
t0 = np.log2(Nsc)

#Calculate time for Brownian oscilations
t = msc*80

#Calculate x for Brownian oscilations plot
x_b = -2*(D*t/Rn)**(-2) + np.log(D*t/Rn)

#Calculate x for large Levy jumps
x_l = np.log(D*(t0+msc*80)/Rn)

#Prepare scatter plot for Brownian oscilations
#p1 = plt.scatter(x_b, y, label = 'Data')
plt.scatter(x_b, y, label = 'Data')

#Linear fit for Brownian oscilations
m, b = np.polyfit(x_b, y, 1)

p2 = plt.plot([np.min(x_b),np.max(x_b)], [m*np.min(x_b),m*np.max(x_b)]+b,
color='r', 
label = 'Linear fit',
linestyle='dashed')
print('Slope for brownian ooscilatiosn m =',+m)

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
plt.scatter(x_l, y, label = 'Data')

#Linear fit for large Levy jumps
m = 1
n = np.mean(y - m*x_l)
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
#plt.show()

