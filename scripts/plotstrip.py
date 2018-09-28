'''
Created on May 25, 2018

@author: aribesme
'''

import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import sys
import os.path
from pylab import rcParams
rcParams['figure.figsize'] = 15, 5
def main():
        result = 'Extraelements/'
	path = 'RESULTS/MDcycle/flattop_std_IP15/'
	ip = '_IP15'
	plotlistfh = ['hoff_f.MQXFA.A1L5','hoff_f.MQXFA.A1R5','hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXFA.A1L1','hoff_f.MQXFA.A1R1','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1']
	plotlistbh = ['hoff_b.MQXFA.A1L5','hoff_b.MQXFA.A1R5','hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXFA.A1L1','hoff_b.MQXFA.A1R1','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1']

	plotlistfv = ['voff_f.MQXFA.A1L5','voff_f.MQXFA.A1R5','voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXFA.A1L1','voff_f.MQXFA.A1R1','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1']
	plotlistbv = ['voff_b.MQXFA.A1L5','voff_b.MQXFA.A1R5','voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXFA.A1L1','voff_b.MQXFA.A1R1','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1']

#  	plotlistfh = ['hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1','hoff_f.MQXFA.A3L1','hoff_f.ACFCA.CL1.B1','hoff_f.ACFCA.CR1.B1','hoff_f.ACFCA.CL5.B1','hoff_f.ACFCA.CR5.B1']
#	plotlistbh = ['hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1','hoff_b.MQXFA.A3L1','hoff_b.ACFCA.CL1.B1','hoff_b.ACFCA.CR1.B1','hoff_b.ACFCA.CL5.B1','hoff_b.ACFCA.CR5.B1']

#	plotlistfv = ['voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1','voff_f.MQXFA.A3L1','voff_f.ACFCA.CL1.B1','voff_f.ACFCA.CR1.B1','voff_f.ACFCA.CL5.B1','voff_f.ACFCA.CR5.B1']
#	plotlistbv = ['voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1','voff_b.MQXFA.A3L1','voff_b.ACFCA.CL1.B1','voff_b.ACFCA.CR1.B1','voff_b.ACFCA.CL5.B1','voff_b.ACFCA.CR5.B1']


	
#	plotlistfh = ['hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1','hoff_f.MQXA.3L2','hoff_f.MQXA.3R2','hoff_f.MQXA.3L8','hoff_f.MQXA.3R8']
#	plotlistbh = ['hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1','hoff_b.MQXA.3L2','hoff_b.MQXA.3R2','hoff_b.MQXA.3L8','hoff_b.MQXA.3R8']

#	plotlistfv = ['voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1','voff_f.MQXA.3L2','voff_f.MQXA.3R2','voff_f.MQXA.3L8','voff_f.MQXA.3R8']
#	plotlistbv = ['voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1','voff_b.MQXA.3L2','voff_b.MQXA.3R2','voff_b.MQXA.3L8','voff_b.MQXA.3R8']
#	plotlistfh = ['hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1']
#	plotlistbh = ['hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1']
#	plotlistfv = ['voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1']
#	plotlistbv = ['voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1']
	

#	plotlistfh = ['hoff_f.ACFCA.CL1.B1','hoff_f.ACFCA.CR1.B1','hoff_f.ACFCA.CL5.B1','hoff_f.ACFCA.CR5.B1']
	
#	plotlistbh = ['hoff_b.ACFCA.CL1.B1','hoff_b.ACFCA.CR1.B1','hoff_b.ACFCA.CL5.B1','hoff_b.ACFCA.CR5.B1']
#	plotlistfv = ['voff_f.ACFCA.CL1.B1','voff_f.ACFCA.CR1.B1','voff_f.ACFCA.CL5.B1','voff_f.ACFCA.CR5.B1']
#	plotlistbv = ['voff_b.ACFCA.CL1.B1','voff_b.ACFCA.CR1.B1','voff_b.ACFCA.CL5.B1','voff_b.ACFCA.CR5.B1']	

	for i in range(len(plotlistbv)):
                print(path+plotlistfh[i])
 		
		bf,fh = np.loadtxt(path+plotlistfh[i],unpack=True, skiprows=0, usecols=[0, 1])
		print(path+plotlistbh[i])
		bb,bh = np.loadtxt(path+plotlistbh[i],unpack=True, skiprows=0, usecols=[0, 1])        
		print(path+plotlistfv[i])
                bf,fv = np.loadtxt(path+plotlistfv[i],unpack=True, skiprows=0, usecols=[0, 1])
		print(path+plotlistbv[i])		
		bb,bv = np.loadtxt(path+plotlistbv[i],unpack=True, skiprows=0, usecols=[0, 1])
		fv = abs(fv)
 		fh = abs(fh)
		bh = abs(bh)
		bv = abs(bv)
                #if int(plotlistbv[i][-1]) == 5 and plotlistbv[i][7] == 'M' :
     		#	fv = abs(fv)
  		#	bv = abs(bv)
 		#elif int(plotlistbv[i][-1]) == 1 and plotlistbv[i][7] == 'M':
		#	fh = abs(fh)
 		#	bh = abs(bh)

		fig  = plt.figure()
		#plt.plot(bf, fh,'rp',label='B1 {st}'.format(st=plotlistfh[i]),alpha=0.5);
		plt.plot(bf, fh,'rp',label='B1',alpha=0.5);
		plt.title('Horizontal offset at {st}'.format(st=plotlistfh[i][7:]));
		plt.xlabel('slot id',fontsize = 18);
		plt.ylabel('Absolute value of the orbit offset (um)',fontsize = 18);
		plt.grid(True);
                plt.xticks(fontsize=15, rotation=0);
		plt.yticks(fontsize=15, rotation=0);
		plt.legend();
		#fig.savefig(result+'{st}'.format(st=plotlistfh[i])+ip+ '_b1.png')
		#plt.show()
		#fig = plt.figure()
                #plt.plot(bb,bh,'bp',label='B2 {st}'.format(st=plotlistfh[i]),alpha=0.5);
                plt.plot(bb,bh,'bp',label='B2',alpha=0.5);
		plt.title('Horizontal offset at {st}'.format(st=plotlistfh[i][7:]));
		plt.xlabel('slot id',fontsize = 18);
		plt.ylabel('Absolute value of the orbit offset (um)',fontsize = 18);
		plt.grid(True);
                plt.xticks(fontsize=15, rotation=0);
		plt.yticks(fontsize=15, rotation=0);
		plt.legend();
		fig.savefig(result+'{st}'.format(st=plotlistfh[i])+ip+ '_absval.png')
		plt.show()
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
		
		#plt.show()
		fig  = plt.figure()
		plt.plot(bf, fv,'rp',label='B1',alpha=0.5)
		plt.title('Vertical offset at {st}'.format(st=plotlistfv[i][7:]))
		plt.xlabel('slot id',fontsize = 18)
		plt.ylabel('Absolute value of the orbit offset (um)',fontsize = 18)
		plt.grid(True);
                plt.xticks(fontsize=15, rotation=0);
		plt.yticks(fontsize=15, rotation=0);
		plt.legend();
		#fig.savefig(result+'{st}'.format(st=plotlistfv[i])+ ip+'_b1.png')
		#plt.show()
                #fig = plt.figure()
                plt.plot(bb,bv,'bp',label='B2',alpha=0.5)
		plt.title('Vertical offset at {st}'.format(st=plotlistfv[i][7:]))
		plt.xlabel('slot id',fontsize = 18)
		plt.ylabel('Absolute value of the orbit offset (um)',fontsize = 18)
		plt.grid(True);
                plt.xticks(fontsize=15, rotation=0);
		plt.yticks(fontsize=15, rotation=0);
		plt.legend();
			#plt.legend(["B1"]);
			#plt.xlim([0.0, 1000.0]);
			#plt.show()
		plt.show() 
		fig.savefig(result+'{st}'.format(st=plotlistfv[i])+ip+ '_absval.png')
		

if __name__ == '__main__':
	main()












