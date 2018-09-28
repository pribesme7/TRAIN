'''
Created on Jan 17, 2018

@author: aribesme
'''

import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import sys
import os.path


def main(arg1):
	version = arg1
	resultFolder = 'MDcycleplots/'
	version = str(version)
	version = version.replace('/',' ').split()
	print('version',version)
	name = version[-1]
	arg2 = 'RESULTS/MDcycle/'
	for i in [1,2,5,8]:
		#print('i',i,arg2 + arg1 + '/hoff_f.IP{st}'.format(st = i))
		if os.path.isfile(arg2 + arg1 + '/hoff_f.IP{st}'.format(st = i)):
			print('Found files for IP{st}'.format(st = i))
			slot1, IPhoffB1 = np.loadtxt(arg2 + arg1 + '/hoff_f.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			slot11, IPvoffB1 = np.loadtxt(arg2 + arg1 + '/voff_f.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			slot2, IPhoffB2 = np.loadtxt(arg2 + arg1 + '/hoff_b.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			slot21, IPvoffB2 = np.loadtxt(arg2 + arg1 + '/voff_b.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			slot1, IPhoffIPsig = np.loadtxt(arg2 + arg1 + '/hsep_sig.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			slot2, IPvoffIPsig = np.loadtxt(arg2 + arg1 + '/vsep_sig.IP{st}'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1]);
			
			fig  = plt.figure()
			plt.plot(slot11, IPhoffB1,'r.')
			plt.title('Horizontal offset at IP{st}, B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('um',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'IP{st}hoffB1'.format(st=i)+ name+ '.png');
			
			fig  = plt.figure()
			plt.plot(slot11, IPvoffB1,'r.')
			plt.title('Vertical offset at IP{st},B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('um',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'IP{st}voffB1'.format(st=i)+ name+ '.png');
		
			fig  = plt.figure()
			plt.plot(slot21, IPhoffB2,'r.')
			plt.title('Horizontal offset at IP{st}, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('um',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B2"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'IP{st}hoffB2'.format(st=i)+ name+ '.png');
	
			fig  = plt.figure()
			plt.plot(slot21, IPvoffB2,'r.')
			plt.title('Vertical offset at IP{st}, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('um',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B2"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'IP{st}voffB2'.format(st=i)+ name+ '.png');
	
			fig  = plt.figure()
			plt.title('Separation at IP{st} '.format(st=i),fontsize = 18);
			plt.plot(slot1, IPhoffIPsig, 'k.');
			plt.plot(slot2, IPvoffIPsig, 'r.');
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('sig',fontsize = 14);
			plt.legend(["Horizontal", "Vertical"]);
			plt.grid(True);
			#plt.show()
			fig.savefig(resultFolder + 'sepIP{st}'.format(st=i)+ name+ '.png');



	
	
if __name__ == "__main__":
    main(sys.argv[1])	






