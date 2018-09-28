# Auxiliar program for writing tables
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import xticks
from pylab import rcParams

rcParams['figure.figsize'] = 20, 5
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16

def read_opt(name):
        path = name + '/fort.fort.ps.hoff_fIP5'
	with open (path,'r') as file:
                cont = 0
    		for line in file:
                        if cont == 0:
        			line = line.split()
       				h1 = round(float(line[-1]),2)
				cont = cont + 1

        path = name + '/fort.fort.ps.voff_fIP5'
	with open (path,'r') as file:
                cont = 0
    		for line in file:
                        if cont == 0:
        			line = line.split()
       				v1 = round(float(line[-1]),2)
				cont = cont + 1


        path = name + '/fort.fort.ps.hoff_bIP5'
	cont = 0
	with open (path,'r') as file:

    		for line in file:
                        if cont == 0:
        			line = line.split()
       				h2 = round(float(line[-1]),2)
				cont = cont + 1

        path = name + '/fort.fort.ps.voff_bIP5'
	cont = 0
	with open (path,'r') as file:
    		for line in file:
			if cont == 0:
	        		line = line.split()
       				v2 = round(float(line[-1]),2)
        			cont = cont + 1 
	return h1,v1,h2,v2	


def read_tune_chrom(name):
        path = name + '/fort.tuneschromashift'
        values = []
	with open (path,'r') as file:
    		for line in file:
        		line = line.split()
        		values.append(float(line[-1]))

        for i in range(len(values)):
		if i == 0 or i == 1 or i == 4 or i ==5:
			values[i] = round(values[i]/10**(-3),3)
		else:
			values[i] = round(values[i],3)
	
	th1 = values[0]
	tv1 = values[1]
        ch1 = values[2]
        cv1 = values[3]
	
	th2 = values[4]
	tv2 = values[5]
        ch2 = values[6]
        cv2 = values[7]
	return th1,tv1,ch1,cv1,th2,tv2,ch2,cv2

def average_chroma(name,skip):
	pathf = name + '/tune_f.list'
	pathb = name + '/tune_b.list'
        chf,cvf = np.loadtxt(pathf,unpack=True,usecols=[4,5],skiprows=skip)
	chb,cvb = np.loadtxt(pathb,unpack=True,usecols=[4,5],skiprows=skip)
	
	avchf = sum(chf)/float(len(chf))
	avcvf = sum(cvf)/float(len(cvf))
	avchb = sum(chb)/float(len(chb))
	avcvb = sum(cvb)/float(len(cvb))
	
	return avchf,avcvf,avchb,avcvb  

skip = [24,12,12]
fscheme = ['std_','B_','8_']
ip_list = ['IP1258']
list_names = ['injection_','rampsqueeze_','presqueeze_','collision_41_','collision_64_','stable_b_41_','stable_b_64_','stable_e_15_']
#folders = ['MDcyclewithoutoct','MDcyclewithoct','MDcycle_ondispon']
folders = ['MDcycle_new_no_oct_no_disp','MDcycle_new_oct_no_disp','MDcycle_new_no_oct_disp','MDcycle_new_oct_disp']
forwardfile = open('table_chroma_forward.dat','w')
backwardfile = open('table_chroma_backward.dat','w')

th1,tv1,ch1,cv1,th2,tv2,ch2,cv2, th1all,tv1all,ch1all,cv1all,th2all,tv2all,ch2all,cv2all = (np.zeros(len(folders)) for i in range(16))
avchf,avcvf,avchb,avcvb = (np.zeros((len(folders),len(list_names),len(fscheme))) for i in range(4))
print(len(avchf))
for k in list_names:
	for l in fscheme:
		for m in range(len(folders)):
			for i in ip_list:
				name = 'RESULTS/{fld}/'.format(fld=folders[m]) + k + l + i
        	                print(name)
        	                if i == 'IP15':
					
					th1[m],tv1[m],ch1[m],cv1[m],th2[m],tv2[m],ch2[m],cv2[m] = read_tune_chrom(name)
				else:
					
					th1all[m],tv1all[m],ch1all[m],cv1all[m],th2all[m],tv2all[m],ch2all[m],cv2all[m] = read_tune_chrom(name)
					
					avchf[m][list_names.index(k)][fscheme.index(l)],avcvf[m][list_names.index(k)][fscheme.index(l)],avchb[m][list_names.index(k)][fscheme.index(l)],avcvb[m][list_names.index(k)][fscheme.index(l)] = average_chroma(name,skip[int(fscheme.index(l))])	
					print(m,k,list_names.index(k),fscheme.index(l),avchf[m][list_names.index(k)][fscheme.index(l)])				
					
		if l == 'std_':
        		fl = 'Standard'
        	elif l == 'B_':
        		fl = 'BCMS'
		else:
			fl = '8b4e'  
		#output1 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12} &{n13}&{n14} &{n15} &{n16} \n'.format(fs = fl,n1 = ch1[0],n2=cv1[0],n3 = ch1all[0],n4=cv1all[0],n5 = ch1[1],n6=cv1[1],n7 = ch1all[1],n8=cv1all[1],n9 = ch1[2],n10=cv1[2],n11 = ch1all[2],n12=cv1all[2],n13 = ch1[3],n14=cv1[3],n15 = ch1all[3],n16=cv1all[3])
		#output1 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12} \n'.format(fs = fl,n1 = ch1[0],n2=cv1[0],n3 = ch1all[0],n4=cv1all[0],n5 = ch1[1],n6=cv1[1],n7 = ch1all[1],n8=cv1all[1],n9 = ch1[2],n10=cv1[2],n11 = ch1all[2],n12=cv1all[2])
		output1 = '~{fs} & {n3}& {n4} &{n7}& {n8} &{n11} &{n12}  &{n5} &{n6} \n'.format(fs = fl,n3 = round(ch1all[0],1),n4=round(cv1all[0],1),n7 = round(ch1all[1],1),n8=round(cv1all[1],1),n11 = round(ch1all[2],1),n12=round(cv1all[2],1),n5 = round(ch1all[3],1),n6=round(cv1all[3],1))
		forwardfile.write(output1)
		
		print(output1)
		#output2 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12} \n'.format(fs = fl,n1 = ch2[0],n2=cv2[0],n3 = ch2all[0],n4=cv2all[0],n5 = ch2[1],n6=cv2[1],n7 = ch2all[1],n8=cv2all[1],n9 = ch2[2],n10=cv2[2],n11 = ch2all[2],n12=cv2all[2])
		output2 = '~{fs} & {n3}& {n4} &{n7}& {n8} &{n11} &{n12} &{n5} &{n6} \n'.format(fs = fl,n3 = round(ch2all[0],1),n4=round(cv2all[0],1),n7 = round(ch2all[1],1),n8=round(cv2all[1],1),n11 = round(ch2all[2],1),n12=round(cv2all[2],1),n5 = round(ch2all[3],1),n6=round(cv2all[3],1))
		
#		output2 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12}&{n13}&{n14} &{n15} &{n16} \n'.format(fs = fl,n1 = ch2[0],n2=cv2[0],n3 = ch2all[0],n4=cv2all[0],n5 = ch2[1],n6=cv2[1],n7 = ch2all[1],n8=cv2all[1],n9 = ch2[2],n10=cv2[2],n11 = ch2all[2],n12=cv2all[2],n13 = ch2[3],n14=cv2[3],n15 = ch2all[3],n16=cv2all[3])
		backwardfile.write(output2)
		print(output2)
	

forwardfile.close()
backwardfile.close()
print(avchf[0,0:,0])
print(len(list_names),len(avchf[0][0:][0]))
plt.plot(range(len(list_names)),avchf[0,0:,0],'o',label=r"$ Q $",markersize=12,alpha=0.5)
plt.plot(range(len(list_names)),avchf[1,0:,0],'o',label=r"$ Q_{fd} $",markersize=12,alpha=0.5)
plt.plot(range(len(list_names)),avchf[2,0:,0],'o',label=r"$ Q_{disp} $",markersize=12,alpha=0.5)
plt.plot(range(len(list_names)),avchf[3,0:,0],'o',label=r"$ Q_{disp+fd} $",markersize=12,alpha=0.5)
xticks(range(7), ('Injection', 'Ramp and Squeeze','Presqueeze','Collision (41cm)','Collision (64cm)', 'Stable beam (41cm)','Stable beam (64cm)'),fontsize = 18)
plt.plot(range(2),20*np.ones(2),'k--')
plt.plot(range(2,5),15*np.ones(3),'k--')
plt.plot(range(5,7),5*np.ones(2),'k--')
plt.ylabel('Average chromaticity',fontsize = 18)
plt.legend()
plt.show()
