# Auxiliar program for writing tables
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import xticks
from pylab import rcParams

rcParams['figure.figsize'] = 20, 5
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16

phase = '/presqueeze_std_IP1258/'
skip =24
filef = 'tune_f.list'
fileb = 'tune_b.list'

folders = ['MDcycle_new_no_oct_no_disp','MDcycle_new_oct_no_disp','MDcycle_new_no_oct_disp','MDcycle_new_oct_disp']
labels = [r"$ Q $",r"$ Q_{fd} $",r"$ Q_{disp} $",r"$ Q_{disp and fd} $"]
for i in folders:
	namef = 'RESULTS/'+i+phase+filef
	nameb = 'RESULTS/'+i+phase+fileb
	nf,thf,tvf,chf,cvf = np.loadtxt(namef,unpack=True,usecols=[0,2,3,4,5],skiprows=skip)
	nb,thb,tvb,chb,cvb = np.loadtxt(nameb,unpack=True,usecols=[0,2,3,4,5],skiprows=skip)
	
	plt.plot(nf,thf,'o',label=labels[folders.index(i)],alpha=0.5)
plt.legend()	
plt.xlabel('Bunch id',fontsize=18)
plt.ylabel('Horizontal tune',fontsize=18)
plt.grid(True)
plt.show()
