# Auxiliar program for writing tables

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
			values[i] = round(values[i]/10**(-3),1)
		else:
			values[i] = round(values[i],1)
	
	th1 = values[0]
	tv1 = values[1]
        ch1 = values[2]
        cv1 = values[3]
	
	th2 = values[4]
	tv2 = values[5]
        ch2 = values[6]
        cv2 = values[7]
	return th1,tv1,ch1,cv1,th2,tv2,ch2,cv2
        


fscheme = ['std_','B_','8_']
ip_list = ['IP15','IP1258']
list_names = ['injection_','flattop_','collision_','squeeze_','collision_ultimate_','stable_b_']

forwardfile = open('table_forward.dat','w')
backwardfile = open('table_backward.dat','w')

for k in list_names:
	for l in fscheme:
		for i in ip_list:
			name = 'RESULTS/MDcyclewithoutoctupoles/' + k + l + i
                        print(name)
                        if i == 'IP15':
				h1,v1,h2,v2 = read_opt(name)
				th1,tv1,ch1,cv1,th2,tv2,ch2,cv2 = read_tune_chrom(name)
			else:
				h1all,v1all,h2all,v2all = read_opt(name)
				th1all,tv1all,ch1all,cv1all,th2all,tv2all,ch2all,cv2all = read_tune_chrom(name)
		if l == 'std_':
        		fl = 'Standard'
        	elif l == 'B_':
        		fl = 'BCMS'
		else:
			fl = '8b4e'  

		output1 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12} \n'.format(fs = fl,n1 = h1,n2 = v1,n3 = h1all,n4 = v1all, n5 = th1,n6 = tv1,n7 = th1all,n8 = tv1all,n9 = ch1,n10 = cv1,n11 = ch1all,n12 = cv1all)
		forwardfile.write(output1)
		print(output1)
		output2 = '~{fs} &{n1} &{n2} & {n3}& {n4} &{n5} &{n6} &{n7}& {n8}&{n9} &{n10} &{n11} &{n12} \n'.format(fs = fl,n1 = h2,n2 = v2,n3 = h2all,n4 = v2all, n5 = th2,n6 = tv2,n7 = th2all,n8 = tv2all,n9 = ch2,n10 = cv2,n11 = ch2all,n12 = cv2all)
		backwardfile.write(output2)
	

forwardfile.close()
backwardfile.close()
