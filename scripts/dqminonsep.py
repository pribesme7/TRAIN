from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import  numpy as np
import subprocess
import pickle as pk

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

    subst5 = 'on_sep5              :=  {st} ;'.format(st=stop )
    pattern5 = 'on_sep5              :=  {st} ;'.format(st=start)
    replace(filename, subst5, pattern5)

    f = open(filename, 'r')
    content = f.read()
    print(content)
    f.close()


def analyse_onsep(outfile,outfile2,onseprange,start,step,filename,prc,stg):

    for onsep in onseprange:
   	    
	    subst = 'on_sep5              :=  {st} ;'.format(st = onsep - step)
            pattern = 'on_sep5              :=  {st} ;'.format(st = onsep)

            replace(filename,subst,pattern)
            
            f = open(filename,'r')
	    content = f.read()
	    print('CHECKING CONTENT')
	    print(content)
            f.close()

            onsepm = int(round((onsep - start - step )/step))

            analyse_dqmin(outfile,outfile2,onsep,prc,stg)

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

    return dqmin,xf0,yf0,xb0,yb0

def analyse_dqmin(outfile,outfile2,onsep,prc,stg):

    file = '25ns_zero_intensity.in'

    subprocess.call(["./updateMadFiles.sh","HL_{pr}".format(pr = prc)])
    
    
    dqmin_analytic,xf0,yf0,xb0,yb0 = analytic_approach()
 
    subprocess.call(["./runTrainForFillingScheme.sh","{f}".format(f = file),"noPlot"])
    print('FINISHED WITH TRAIN, WRITING THE RESULTS!!!!')
    filepathf = "RESULTS/{f}/fort.closest_tune_app_f".format(f = file)

    filepathb = "RESULTS/{f}/fort.closest_tune_app_b".format(f = file)


    ff = open(filepathf,'r')
    fb = open(filepathb,'r')

    with open(filepathf,'r') as ff:
	with open(filepathb,'r') as fb:
    		aux, dqminf,betxf,betyf,gammaf,rf11,rf12,rf21,rf22 = np.loadtxt(ff, unpack=True, skiprows=0, usecols=[0, 1,2,3,4,5,6,7,8])
    		aux, dqminb,betxb,betyb,gammab,rb11,rb12,rb21,rb22 = np.loadtxt(fb, unpack=True, skiprows=0, usecols=[0, 1,2,3,4,5,6,7,8])	

    

    		outcontentf = '{d} \t {xf} \t {yf} \t {dqmf} \t {dqmanl} \t {btxf} \t {btyf} \t {gmmf} \t {rf11} \t {rf12} \t {rf21} \t {rf22} \n'.format(d = onsep,xf = xf0, yf = yf0, dqmf = dqminf,dqmanl = dqmin_analytic,btxf = betxf,btyf = betyf,gmmf = gammaf,rf11 = rf11, rf12 = rf12,rf21 = rf21,rf22 = rf22)
    		outcontentb = '{d} \t {xb} \t {yb} \t {dqmb} \t {dqmanl} \t {btxb} \t {btyb} \t {gmmb} \t {rb11} \t {rb12} \t {rb21} \t {rb22} \n'.format(d = onsep,xb = xb0, yb = yb0, dqmb = dqminb,dqmanl = dqmin_analytic,btxb = betxb,btyb = betyb,gmmb = gammab,rb11 = rb11, rb12 = rb12,rb21 = rb21,rb22 = rb22)
                print(outcontentf)
                print(outcontentb)

    		outfile.write(outcontentf)
    		outfile2.write(outcontentb)
    
    

def main():
    init = False
    prc = "diagonal"
    stg = "zero"

    step = 0.5
    start = -20.0 - step
    stop = 20.0 + step
 
    

    filename = "MAD_PART/MDcycle/on_sep.{pr}".format(pr = prc)

    if init == True:
        initialize(filename, start, 20.0)
	
    else:
	onsep = np.arange(start + step,-step,step)
	onsep2 = np.arange(-step,stop,step)


	onsep = list(onsep) + list(onsep2)
	onsep = [ '%.2f' % elem for elem in onsep ]
	onseprange = [ float(elem) for elem in onsep ]
	
	for i in onseprange:
		print(int(round((i - start - step )/step)))
        

	filename2 = 'RESULTS/dqminonsep_f.dat'
        filename3 = 'RESULTS/dqminonsep_b.dat'
	outfile = open(filename2,'w')
        outfile2 = open(filename3,'w')
#        with open(filename2, 'w') as outfile:
#		with open(filename3, 'w') as outfile2:
#      	outfile.write('# Linear coupling forward beam \n')
#        outfile2.write('# Linear coupling backward beam \n')
        analyse_onsep(outfile,outfile2,onseprange,start,step,filename,prc,stg)

        outfile.close()
        outfile2.close()
if __name__ == "__main__":
    main()
