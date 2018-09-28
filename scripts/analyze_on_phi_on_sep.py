from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import  numpy as np
import subprocess
import pickle as pk

def parse_tfs_from_MADX(filename,delimiter=" ",header_identifier="@",categories_startswith="*",descriptor_startswith="$"):
   
    with open(filename, "r") as f:
        lines = f.readlines()
    i = -1
    data = {}
    header = {}
    descriptor_line =[]
    keys = []
    for line in lines: 
        line = filter(None,line.strip('\n').split(delimiter))
         
        if line[0] == header_identifier:
            try:
                header[line[1]] = float(line[3])
            except:
                header[line[1]] = " ".join(line[3:])
        elif line[0] ==categories_startswith:
            for l in line[1:]: 
                data[l] = []
            keys = line[1:]
        elif line[0] == descriptor_startswith:
            descriptor_line = line[1:]
        else: 
            for j,l in enumerate(line):
               try: 
                    data[keys[j]].append(float(l))
               except:
                   data[keys[j]].append(l)
    
    return header, data 


def analytic_approach():
    tf = 'train.optf'
    tb = 'train.optb'
    

    headerf,summf = parse_tfs_from_MADX(tf,delimiter=" ",header_identifier="@",categories_startswith="*",descriptor_startswith="$")
    headerb,summb = parse_tfs_from_MADX(tb,delimiter=" ",header_identifier="@",categories_startswith="*",descriptor_startswith="$")
    



    namef = summf['NAME']
    nameb = summb['NAME']

    xf = summf['X']
    xb = summb['X']

    xf0 = xf[namef.index('"IP5"')]
    xb0 = xb[nameb.index('"IP5"')]

    betxf = summf['BETX']
    betxb = summb['BETX']

    betxf0 = betxf[namef.index('"IP5"')]
    betxb0 = betxb[namef.index('"IP5"')]    

    yf = summf['Y']
    yb = summb['Y']

    yf0 = yf[namef.index('"IP5"')]
    yb0 = yb[nameb.index('"IP5"')]
 
    betyf = summf['BETY']
    betyb = summb['BETY']

    betyf0 = betyf[namef.index('"IP5"')]
    betyb0 = betyb[namef.index('"IP5"')] 
   
    r0 = 1.53469857 * 10 ** (-18)
    en = 2.5 *10**(-6)
    gamma = float(headerf['GAMMA'])
    N = float(headerf['NPART'])
    
    sigmaxf = np.sqrt(en*betxf0/gamma)
    sigmayf = np.sqrt(en*betyf0/gamma)
    sigmaf = np.sqrt(sigmaxf**2 + sigmayf**2)
    
    sigmaxb = np.sqrt(en*betxb0/gamma)
    sigmayb = np.sqrt(en*betyb0/gamma)
    sigmab = np.sqrt(sigmaxb**2 + sigmayb**2) 
    sigma = np.sqrt(sigmaf**2 + sigmab**2)
    d = np.sqrt((xf0-xb0)**2+(yf0-yb0)**2)/sigma
    alpha = abs(np.arctan((yf0-yb0)/(xf0-xb0)))

    dqmin = 2*r0*N/(np.pi*en)*np.sin(alpha)*np.cos(alpha)/d**2
    dqtrial = r0 * N / (2*np.pi * en) * np.sin(alpha) * np.cos(alpha) / d ** 2
    return dqmin,dqtrial


def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)


def initialize(filename,start,stop):

    subst1 = 'phi_IR1              :=   {st} ;'.format(st=stop )
    pattern1 = 'phi_IR1              :=   {st} ;'.format(st=start)
    replace(filename, subst1, pattern1)

    subst5 = 'phi_IR5              :=   {st} ;'.format(st=stop )
    pattern5 = 'phi_IR5              :=   {st} ;'.format(st=start)
    replace(filename, subst5, pattern5)

    f = open(filename, 'r')
    content = f.read()
    print(content)
    f.close()


def analyse_phi(outfile,philist,step,filemodify,ip):

	for phi in philist:
   	    
	    subst = 'phi_IR1              :=   {st} ;'.format(st = phi - step)
            pattern = 'phi_IR1              :=   {st} ;'.format(st = phi    )

            replace(filemodify,subst,pattern)
            
            f = open(filemodify,'r')
	    content = f.read()
	    print('CHECKING CONTENT')
	    print(content)
            f.close()
	    
	    execute_train(outfile,phi,ip) 


def analyse_2phi(outfile,phi_1,phi_5,step,filemodify):

	for phi1 in phi_1:
		subst1 = 'phi_IR1              :=   {st} ;'.format(st = phi1 - step)
            	pattern1 = 'phi_IR1              :=   {st} ;'.format(st = phi1)
		replace(filemodify,subst1,pattern1)
		for phi5 in phi_5[0:-1]:
			print(phi1,phi5)
			subst5 = 'phi_IR5              :=   {st} ;'.format(st = phi5)
            		pattern5 = 'phi_IR5              :=   {st} ;'.format(st = phi5 - step)
			
			replace(filemodify,subst5,pattern5)
            		print(phi1,phi5)
            
            		f = open(filemodify,'r')
	    		content = f.read()
	    		print('CHECKING CONTENT')
	    		print(content)
            		f.close()
			phi2 = phi5 - step
			execute_train2(outfile,phi1,phi2)

	    	subst5 = 'phi_IR5              :=   {st} ;'.format(st = 45)
            	pattern5 = 'phi_IR5              :=   {st} ;'.format(st = 91)
		replace(filemodify,subst5,pattern5)
		
	    		 

def execute_train2(outfile,phi1,phi5):

	collision41 = "HL_collision_41"
	std = "25ns_2760b_2748_2494_2572_288bpi_13inj.in"

	subprocess.call(["./updateMadFiles.sh",collision41])
	subprocess.call(["./runTrainForFillingScheme.sh","{stf}".format(stf=std),"noPlot"])
	
	filepathf = "RESULTS/{stf}/fort.closest_tune_app_f".format(stf = std)
    	filepathb = "RESULTS/{stf}/fort.closest_tune_app_b".format(stf = std)

	filepath1 = "RESULTS/{stf}/fort.fort.ps.hoff_fIP5".format(stf=std)
	filepath2 = "RESULTS/{stf}/fort.fort.ps.hoff_bIP5".format(stf=std)
	filepath3 = "RESULTS/{stf}/fort.fort.ps.voff_fIP5".format(stf=std)
	filepath4 = "RESULTS/{stf}/fort.fort.ps.voff_bIP5".format(stf=std)
	filepath5 = "RESULTS/{stf}/fort.tuneschromashift".format(stf=std)
	
	mtaf = np.loadtxt(filepathf,unpack=True,usecols=[1],skiprows=24)
	mtab = np.loadtxt(filepathb,unpack=True,usecols=[1],skiprows=24)

	hf = np.loadtxt(filepath1,unpack=True,usecols=[-1])
	hb = np.loadtxt(filepath2,unpack=True,usecols=[-1])
	vf = np.loadtxt(filepath3,unpack=True,usecols=[-1])
	vb = np.loadtxt(filepath4,unpack=True,usecols=[-1])
	tc = np.loadtxt(filepath5,unpack=True,usecols=[-1])
	
	dqspread1 = max(mtaf) - min(mtaf)
	dqaverage1 = sum(mtaf)/float(len(mtaf))
	
	dqspread2 = max(mtab) - min(mtab)
	dqaverage2 = sum(mtab)/float(len(mtab))	
	
	dq,dqtrial = analytic_approach()	
	
	output = "{s01} \t {s0} \t {s1} \t {s2} \t {s3} \t {s4} \t{ss1} \t {ss2} \t {s5} \t {s6} \t {s7} \t {s8} \t {s9} \t {s10} \t {s11} \t {s12} \t {s13} \t {s14} \t {s15} \t {s16} \n".format(s01=phi1,s0 = phi5,s1=dqaverage1,s2=dqspread1,s3=dqaverage2,s4=dqspread2,ss1=dq,ss2=dqtrial,s5=hf[0],s6=hb[0],s7=vf[0],s8=vb[0],s9=hf[1],s10=hb[1],s11=vf[1],s12=vb[1],s13=tc[0],s14=tc[4],s15=tc[1],s16=tc[5])

	outfile.write(output)
	print(output)


def compute_angle(bets,ip):
	marker1 = []
 	marker2 = []
	with open('train.optf','r') as file1:
		for line in file1:
			line = line.strip().split()
			if line[0] != '@' and line[0] != '*' and line[0]!= '$':
				marker1.append(line[0])
	with open('train.optb','r') as file2:
		for line in file2:
			line = line.strip().split()
			if line[0] != '@' and line[0] != '*' and line[0]!= '$':
				marker2.append(line[0])

	s1,x1,y1 = np.loadtxt('train.optf',unpack=True,usecols=[1,2,6],skiprows=47)
	s2,x2,y2 = np.loadtxt('train.optb',unpack=True,usecols=[1,2,6],skiprows=47)
	s1 = list(s1)
	s2 = list(s2)
	
	sigma = np.sqrt(bets*2.6*10**(-6)/7.460522527611191435e+03)
	d1 = np.sqrt((x1[marker1.index('"MKIP1PL1"')]-x2[marker2.index('"MKIP1PL1"')])**2+(y1[marker1.index('"MKIP1PL1"')]-y2[marker2.index('"MKIP1PL1"')])**2)/sigma
	alpha1 = np.arctan((y1[marker1.index('"MKIP1PL1"')]-y2[marker2.index('"MKIP1PL1"')])/x1[marker1.index('"MKIP1PL1"')]-x2[marker2.index('"MKIP1PL1"')])*180./np.pi

	d5 = np.sqrt((x1[marker1.index('"MKIP5PL1"')]-x2[marker2.index('"MKIP5PL1"')])**2+(y1[marker1.index('"MKIP5PL1"')]-y2[marker2.index('"MKIP5PL1"')])**2)/sigma
	alpha5 = np.arctan((y1[marker1.index('"MKIP5PL1"')]-y2[marker2.index('"MKIP5PL1"')])/x1[marker1.index('"MKIP5PL1"')]-x2[marker2.index('"MKIP5PL1"')])*180./np.pi


	if ip == 5:
		return alpha5
	else:
		return alpha1

def execute_train(outfile,phi,ip):

	collision41 = "HL_collision_41"
	std = "25ns_2760b_2748_2494_2572_288bpi_13inj.in"
        bets = 0.41
	subprocess.call(["./updateMadFiles.sh",collision41])
	subprocess.call(["./runTrainForFillingScheme.sh","{stf}".format(stf=std),"noPlot"])
	
	filepathf = "RESULTS/{stf}/fort.closest_tune_app_f".format(stf = std)
    	filepathb = "RESULTS/{stf}/fort.closest_tune_app_b".format(stf = std)

	filepath1 = "RESULTS/{stf}/fort.fort.ps.hoff_fIP5".format(stf=std)
	filepath2 = "RESULTS/{stf}/fort.fort.ps.hoff_bIP5".format(stf=std)
	filepath3 = "RESULTS/{stf}/fort.fort.ps.voff_fIP5".format(stf=std)
	filepath4 = "RESULTS/{stf}/fort.fort.ps.voff_bIP5".format(stf=std)
	filepath5 = "RESULTS/{stf}/fort.tuneschromashift".format(stf=std)
        
        
	
	mtaf = np.loadtxt(filepathf,unpack=True,usecols=[1],skiprows=24)
	mtab = np.loadtxt(filepathb,unpack=True,usecols=[1],skiprows=24)

	hf = np.loadtxt(filepath1,unpack=True,usecols=[-1])
	hb = np.loadtxt(filepath2,unpack=True,usecols=[-1])
	vf = np.loadtxt(filepath3,unpack=True,usecols=[-1])
	vb = np.loadtxt(filepath4,unpack=True,usecols=[-1])
	tc = np.loadtxt(filepath5,unpack=True,usecols=[-1])
	
	dqspread1 = max(mtaf) - min(mtaf)
	dqaverage1 = sum(mtaf)/float(len(mtaf))
	
	dqspread2 = max(mtab) - min(mtab)
	dqaverage2 = sum(mtab)/float(len(mtab))	
	
	dqanltcav,dqanltcspread =  analytic_cmin(std,ip)

        phi2 = compute_angle(bets,ip)

	output = " {s0} \t {sm2} \t {sm1} \t{s1} \t {s2} \t {s3} \t {s4} \t {s5} \t {s6} \t {s7} \t {s8} \t {s9} \t {s10} \t {s11} \t {s12} \t {s13} \t {s14} \t {s15} \t {s16} \t {s17} \n".format(s0 = phi,sm2 = dqanltcav, sm1 = dqanltcspread, s1=dqaverage1,s2=dqspread1,s3=dqaverage2,s4=dqspread2,s5=hf[0],s6=hb[0],s7=vf[0],s8=vb[0],s9=hf[1],s10=hb[1],s11=vf[1],s12=vb[1],s13=tc[0],s14=tc[4],s15=tc[1],s16=tc[5],s17 = phi2)

	outfile.write(output)
	print(output)

def execute_train_onsep(outfile,phi):

	collision41 = "HL_collision_41"
	std = "25ns_2760b_2748_2494_2572_288bpi_13inj.in"

	subprocess.call(["./updateMadFiles.sh",collision41])
	subprocess.call(["./runTrainForFillingScheme.sh","{stf}".format(stf=std),"noPlot"])
	filename1 = "RESULTS/{stf}/tune_f.list".format(stf=std)
	filename2 = "RESULTS/{stf}/tune_b.list".format(stf=std)
	filepath5 = "RESULTS/{stf}/fort.tuneschromashift".format(stf=std)
	tc = np.loadtxt(filepath5,unpack=True,usecols=[-1])
	bnch1,bcurr1,t1h,t1v,c1h,c1v = np.loadtxt(filename1,unpack=True,usecols=[0,1,2,3,4,5],skiprows = 24)
	bnch2,bcurr2,t2h,t2v,c2h,c2v = np.loadtxt(filename2,unpack=True,usecols=[0,1,2,3,4,5],skiprows = 24)
	
	t1hav = sum(t1h)/float(len(t1h))
	t1vav = sum(t1v)/float(len(t1h))
	t2hav = sum(t2h)/float(len(t1h))
	t2vav = sum(t2v)/float(len(t1h))

	c1hav = sum(c1h)/float(len(c1h))
	c1vav = sum(c1v)/float(len(c1h))
	c2hav = sum(c2h)/float(len(c1h))
	c2vav = sum(c2v)/float(len(c1h))

	output = " {s0} \t {sm2} \t {sm1} \t{s1} \t {s2} \t {s3} \t {s4} \t {s5} \t {s6}\t {s7} \t {s8}\t {s9} \t {s10}\t {s11} \t {s12}\t {s13} \t {s14} \n ".format(s0 = phi,sm2 = t1hav, sm1 = tc[0],s1 = t1vav, s2 = tc[1],s3 = c1hav, s4 = tc[2],s5 = c1vav, s6 = tc[3],s7 = t2hav, s8 = tc[4],s9 = t2vav, s10 = tc[5],s11 = c2hav, s12 = tc[6],s13 = c2vav, s14 = tc[7])	
	outfile.write(output)
	print(output)

def comp_dqmin(npart,gamma,betxf0,betyf0, betxb0,betyb0,xf0,xb0,yf0,yb0):
    r0 = 1.53469857 * 10 ** (-18)
    en = 2.5 * 10 ** (-6)


    sigmaxf = np.sqrt(en * betxf0 / gamma)
    sigmayf = np.sqrt(en * betyf0 / gamma)
    sigmaf = np.sqrt(sigmaxf ** 2 + sigmayf ** 2)

    sigmaxb = np.sqrt(en * betxb0 / gamma)
    sigmayb = np.sqrt(en * betyb0 / gamma)
    sigmab = np.sqrt(sigmaxb ** 2 + sigmayb ** 2)
    sigma = np.sqrt(sigmaf ** 2 + sigmab ** 2)
    print('SIGMA',sigma*1000)
    dx =  xf0 - xb0
    dy = yf0 - yb0
    d = np.sqrt((dx) ** 2 + (dy) ** 2)/sigma
    
    alpha = abs(np.arctan(dy / dx))
    alpha = np.pi/2-alpha
    #if dy < 0 :
#	if dx < 0:
#		alpha = alpha + np.pi
#	else:
#		alpha = -alpha
#    else:
#	if dx < 0:
#		alpha = np.pi - alpha
#    if alpha < 0:
#	alpha = alpha + 2*np.pi

    print('d,alpha',d,dy,dx,alpha*180/np.pi)
    dqmin =  2*r0 * npart / (np.pi * en) * np.sin(alpha) * np.cos(alpha) / d ** 2
    dqtrial = r0 * npart / (2*np.pi * en) * np.sin(alpha) * np.cos(alpha) / d ** 2
    print(dqmin,sigma*1000*20)
    return dqtrial	

def analytic_cmin(std,ip):
	filepath1 = "RESULTS/{stf}/hoff_f.{sip}".format(stf=std,sip=ip)	
	filepath2 = "RESULTS/{stf}/hoff_b.{sip}".format(stf=std,sip=ip)
	filepath3 = "RESULTS/{stf}/voff_f.{sip}".format(stf=std,sip=ip)
	filepath4 = "RESULTS/{stf}/voff_b.{sip}".format(stf=std,sip=ip)
	
	hf = np.loadtxt(filepath1,unpack=True,usecols=[1],skiprows=24)
	hb = np.loadtxt(filepath2,unpack=True,usecols=[1],skiprows=24)
	vf = np.loadtxt(filepath3,unpack=True,usecols=[1],skiprows=24)
	vb = np.loadtxt(filepath4,unpack=True,usecols=[1],skiprows=24)

	betas = 0.41
	gamma = 7.460522527611191435e+03
	npart = 2.3*10**11
	cminus = []
	tm6 = 10**(-6)
	for i in range(len(hf)):
		cminus.append(comp_dqmin(npart,gamma,betas,betas,betas,betas,hf[i]*tm6,hb[i]*tm6,vf[i]*tm6,vb[i]*tm6))
	
	cminspread = max(cminus) - min(cminus)	
	cminav = sum(cminus)/float(len(cminus))

	return cminav,cminspread

def analyse_onsep(outfile,philist,step,filemodify,ip):
	for phi in philist:
   	    
	    subst = 'on_sep1              :=  %.2f ;'% (phi)
            pattern = 'on_sep1              :=  %.2f ;'% (phi- step ) 

            replace(filemodify,subst,pattern)
            
            f = open(filemodify,'r')
	    content = f.read()
	    print('CHECKING CONTENT')
	    print(content)
            f.close()
	    
	    #execute_train_onsep(outfile,phi)	
	    execute_train(outfile,phi,ip)
def main():
	start = 0
	stop = 90
	step = 0.1
	init = False
	filemodify = 'MAD_PART/MDcycle/on_phi.collision41'
        filemodify2 = 'MAD_PART/MDcycle/on_sep.collision'
	prc = "collision41"
        ip = 'IP1'
	fillsch = "std"
	if init == True:
        	initialize(filemodify, start, 20.0)
	
	else:	
		#phi_5 = range(45,stop+2)
		#phi_1 = range(46) 
		#print(phi_5[::-1])
		#phi_5 = phi_5[::-1]
		phi = np.arange(91,45,-1,dtype='int')
		#print(phi)
		outfilename = "{pr}_{fs}_sepscan_lc_ip15_{sip}.dat".format(pr=prc,fs=fillsch,sip=ip)
		outfile = open(outfilename,'w')
           	phi = np.arange(0,46,1,dtype = 'int')
		onsep = np.arange(0.55,-0.55,-0.1)
                print(onsep)
		#phi = np.arange(45,91,1,dtype = 'int')
		analyse_onsep(outfile,onsep,step,filemodify2,ip)	
		#analyse_2phi(outfile,phi_1,phi_5,step,filemodify)
		outfile.close()
		
			

if __name__ == "__main__":
    main()
