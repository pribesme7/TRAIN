import numpy as np
import matplotlib.pyplot as plt
import pickle as p
cycle = ['injection','flattop','collision','squeeze','stable_b']
cycle_nom = ['injection', 'rampsqueeze','collision_64','stable_b_64','stable_e_15']
cycle_ult = ['injection', 'rampsqueeze','presqueeze','collision_41','stable_b_41','stable_e_15']
#fillsch = ['std','B','8']
#strt = [24,12,12]
fillsch = ['8','B','std']
strt = [12,12,24]
ip = 'IP1258/'
path = 'RESULTS/MDcycle/'
#plotlistfh = ['hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1','hoff_f.MQXFA.A3L1','hoff_f.ACFCA.CL1.B1','hoff_f.ACFCA.CR1.B1','hoff_f.ACFCA.CL5.B1','hoff_f.ACFCA.CR5.B1']
#plotlistbh = ['hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1','hoff_b.MQXFA.A3L1','hoff_b.ACFCA.CL1.B1','hoff_b.ACFCA.CR1.B1','hoff_b.ACFCA.CL5.B1','hoff_b.ACFCA.CR5.B1']

#plotlistfv = ['voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1','voff_f.MQXFA.A3L1','voff_f.ACFCA.CL1.B1','voff_f.ACFCA.CR1.B1','voff_f.ACFCA.CL5.B1','voff_f.ACFCA.CR5.B1']
#plotlistbv = ['voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1','voff_b.MQXFA.A3L1','voff_b.ACFCA.CL1.B1','voff_b.ACFCA.CR1.B1','voff_b.ACFCA.CL5.B1','voff_b.ACFCA.CR5.B1']
plotlistfh = ['hoff_f.MQXA.3L8','hoff_f.MQXA.3R8','hoff_f.MQXFA.A3L5','hoff_f.MQXFA.A3R5','hoff_f.MQXA.3L2','hoff_f.MQXA.3R2','hoff_f.MQXFA.A3L1','hoff_f.MQXFA.A3R1','hoff_f.ACFCA.CL1.B1','hoff_f.ACFCA.CR1.B1','hoff_f.ACFCA.CL5.B1','hoff_f.ACFCA.CR5.B1']
plotlistbh = ['hoff_b.MQXA.3L8','hoff_b.MQXA.3R8','hoff_b.MQXFA.A3L5','hoff_b.MQXFA.A3R5','hoff_b.MQXA.3L2','hoff_b.MQXA.3R2','hoff_b.MQXFA.A3L1','hoff_b.MQXFA.A3R1','hoff_b.ACFCA.CL1.B1','hoff_b.ACFCA.CR1.B1','hoff_b.ACFCA.CL5.B1','hoff_b.ACFCA.CR5.B1']

plotlistfv = ['voff_f.MQXA.3L8','voff_f.MQXA.3R8','voff_f.MQXFA.A3L5','voff_f.MQXFA.A3R5','voff_f.MQXA.3L2','voff_f.MQXA.3R2','voff_f.MQXFA.A3L1','voff_f.MQXFA.A3R1','voff_f.ACFCA.CL1.B1','voff_f.ACFCA.CR1.B1','voff_f.ACFCA.CL5.B1','voff_f.ACFCA.CR5.B1']
plotlistbv = ['voff_b.MQXA.3L8','voff_b.MQXA.3R8','voff_b.MQXFA.A3L5','voff_b.MQXFA.A3R5','voff_b.MQXA.3L2','voff_b.MQXA.3R2','voff_b.MQXFA.A3L1','voff_b.MQXFA.A3R1','voff_b.ACFCA.CL1.B1','voff_b.ACFCA.CR1.B1','voff_b.ACFCA.CL5.B1','voff_b.ACFCA.CR5.B1']


pdic = {}
print(pdic)
tempdic = {}
for i in cycle:
	for j in fillsch:
		tempdic[j] = {'h':1,'v':2}
     		
	pdic['b1'] = tempdic

print(pdic)
pdict = {'B1':'','B2':''}
tmpdict1 = {'std':{},'B':{},'8':{}}
tmpdict2 = {'std':{},'B':{},'8':{}}


for k in range(len(plotlistfv)):
	
	for j in range(len(fillsch)):
        
		tmplistfh = []
		tmplistbh = []
		tmplistfv = []
		tmplistbv = []
		
		for i in cycle:
			
			foldname = path + '{s1}_{s2}_'.format(s1=i,s2=fillsch[j])+ip

			bf,fh = np.loadtxt(foldname+plotlistfh[k],unpack=True, skiprows=0, usecols=[0, 1])
			bb,bh = np.loadtxt(foldname+plotlistbh[k],unpack=True, skiprows=0, usecols=[0, 1])        
	                bf,fv = np.loadtxt(foldname+plotlistfv[k],unpack=True, skiprows=0, usecols=[0, 1])
			bb,bv = np.loadtxt(foldname+plotlistbv[k],unpack=True, skiprows=0, usecols=[0, 1])
			
			
			pickfh = max(fh[strt[j]:]) - min(fh[strt[j]:])
 			pickbh = max(bh[strt[j]:]) - min(bh[strt[j]:])
			pickfv = max(fv[strt[j]:]) - min(fv[strt[j]:])
			pickbv = max(bv[strt[j]:]) - min(bv[strt[j]:])
			
			tmplistfh.append(pickfh)
			tmplistbh.append(pickbh)
			tmplistfv.append(pickfv)
			tmplistbv.append(pickbv)
 		tmpdict1[fillsch[j]][plotlistfv[k][7:]] = {'h':tmplistfh,'v':tmplistfv}
 		tmpdict2[fillsch[j]][plotlistfv[k][7:]] = {'h':tmplistbh,'v':tmplistbv}
        
 		
pdict['B1'] = tmpdict1
pdict['B2'] = tmpdict2 		


p.dump( pdict, open( "save_1258.p", "wb" ) )

