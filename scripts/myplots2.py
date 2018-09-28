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
	resultFolder = 'MDcyclenote/'
	version = str(version)
	version = version.replace('/',' ').split()
	print('version',version)
	name = version[-1]
	arg2 = 'RESULTS/MDcycle/'
	for i in ['f','b']:
		#print('i',i,arg2 + arg1 + '/hoff_f.IP{st}'.format(st = i))
		if os.path.isfile(arg2 + arg1 + '/tune_{st}.list'.format(st = i)):
			print('Found files for IP{st}'.format(st = i))
			
			slot1,bcurr1,qx1,qy1,qxp1,qyp1 = np.loadtxt(arg2 + arg1+ '/tune_{st}.list'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1,2,3,4,5]);
			slot2,bcurr2,qx2,qy2,qxp2,qyp2 = np.loadtxt(arg2 + arg1+ '/tune_{st}.list'.format(st = i), unpack=True, skiprows=0, usecols=[0, 1,2,3,4,5]);

			fig  = plt.figure()
			plt.plot(slot1, qx1,'r.')
			plt.title('Horizontal Tune, B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qx1'+ name+ '.png');
			#plt.show()


			fig  = plt.figure()
			plt.plot(slot2, qx2,'r.')
			plt.title('Horizontal Tune, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qx2'+ name+ '.png');
			#plt.show()

			fig  = plt.figure()
			plt.plot(slot1, qy1,'r.')
			plt.title('Vertical Tune, B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qy1'.format(st=i)+ name+ '.png');
			#plt.show()


			fig  = plt.figure()
			plt.plot(slot2, qy2,'r.')
			plt.title('Vertical Tune, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qy2'.format(st=i)+ name+ '.png');
			#plt.show()





			fig  = plt.figure()
			plt.plot(slot1, qxp1,'r.')
			plt.title('Horizontal Chromaticity, B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qxp1'.format(st=i)+ name+ '.png');
			#plt.show()


			fig  = plt.figure()
			plt.plot(slot2, qxp2,'r.')
			plt.title('Horizontal Chromaticity, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qxp2'.format(st=i)+ name+ '.png');
			#plt.show()

			fig  = plt.figure()
			plt.plot(slot1, qyp1,'r.')
			plt.title('Vertical Chromaticity, B1'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qyp1'.format(st=i)+ name+ '.png');
			#plt.show()


			fig  = plt.figure()
			plt.plot(slot2, qyp2,'r.')
			plt.title('Vertical Chromaticity, B2'.format(st=i),fontsize = 18);
			plt.xlabel('slot id',fontsize = 14);
			plt.ylabel('',fontsize = 14);
			plt.grid(True);
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
			fig.savefig(resultFolder + 'qyp2'+ name+ '.png');
			#plt.show()




	
	
if __name__ == "__main__":
    main(sys.argv[1])	

