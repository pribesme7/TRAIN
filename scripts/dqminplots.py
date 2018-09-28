import matplotlib.pyplot as plt
import numpy as np

f = open('RESULTS/dqminonsep_f.dat','r')
b = open('RESULTS/dqminonsep_b.dat','r')
d,x0f,y0f,dqminf,dqant = np.loadtxt(f,unpack=True, skiprows=0, usecols=[0,1,2,3,4])
d,x0b,y0b,dqminb,dqant = np.loadtxt(b,unpack=True, skiprows=0, usecols=[0,1,2,3,4])

plt.plot(d,dqminf,'bo',label = 'Minimum tune approach')
plt.plot(d,dqant,'bo',label = 'Analitic approach')
#plt.legend()
plt.show()
