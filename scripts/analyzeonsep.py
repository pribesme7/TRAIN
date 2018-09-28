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



#filename = "on_sep.collision"
#with open(filename,'r+') as file:
#    for line in file:
#        if line[0] != '!':
#            if int(line[6]) == 1:

#pattern = 'on_sep1              :=  0.0 ;'
#subst = 'on_sep1              :=  {st} ;'.format(st = 1.0)
#replace(filename,subst,pattern)

def initialize(filename,start,stop):

    subst1 = 'on_sep1              :=  {st} ;'.format(st=stop )
    pattern1 = 'on_sep1              :=  {st} ;'.format(st=start)
    replace(filename, subst1, pattern1)

    subst5 = 'on_sep5              :=  {st} ;'.format(st=stop )
    pattern5 = 'on_sep5              :=  {st} ;'.format(st=start)
    replace(filename, subst5, pattern5)

    f = open(filename, 'r')
    content = f.read()
    print(content)
    f.close()



def execute_train(n,outfile,dicip1,dicip5,tiph,tipv,ciph,cipv,onsep1,onsep5,process,stg):
    #on1 = "{onsep1}".format(onsep1= onsep1)
    #on5 ="{onsep5}".format(onsep5=onsep5)
    #subprocess.call(["./prova.sh", str(on1),str(on5)])
    #filename = "on_sep.collision"
    #f = open(filename, 'r')
    #content = f.read()
    #print(content)
    #f.close()
    if stg == 'std':
        file = "25ns_2760b_2748_2494_2572_288bpi_13inj.in"
        nhead = 108
        ntail = 416
        nmid = 483

    elif stg == 'B':
        file = "25ns_2748b_2736_2258_2374_288bpi_12inj.in"
        nhead = 427
        ntail = 395
        nmid = 448

    else:
        file = "8b4e_1972b_1967_1178_1886_224bpi_12inj.in"
        nhead = 420
        ntail = 388
        nmid = 420

    

    subprocess.call(["./updateMadFiles.sh","HL_{pr}".format(pr = process)])

    subprocess.call(["./runTrainForFillingScheme.sh","{f}".format(f = file),"noPlot"])

    subprocess.call(["./subfolders.sh","onsep_{pr}_{st}_IP{num}".format(pr = process,st = stg,num = n)])


# Vertical orbit offset for IP1 
    if int(str(n)[0]) == 1:
    	path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/voff_f.IP1".format(pr = process,st = stg,num=n)
    	f = open(path,'r')

    	col1, col2 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1])
    	col1 = list(col1)
    	col2 = list(col2)
    	dicip1["f"]["head"][onsep1,onsep5] = float(col2[col1.index(nhead)])
    	dicip1["f"]["medium"][onsep1, onsep5] = float(col2[col1.index(nmid)])
    	dicip1["f"]["tail"][onsep1, onsep5] = float(col2[col1.index(ntail)])
    	print(onsep1,onsep5,col1.index(nhead),col1.index(nmid),col1.index(ntail),float(col2[col1.index(nhead)]),float(col2[col1.index(nmid)]),float(col2[col1.index(ntail)]))
    	f.close()

    

    	path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/voff_b.IP1".format(pr = process,st = stg,num=n)
    	f = open(path, 'r')

    	col3, col4 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1])
    	col3 = list(col3)
    	col4 = list(col4)
    	dicip1["b"]["head"][onsep1, onsep5] = float(col4[col3.index(nhead)])
    	dicip1["b"]["medium"][onsep1, onsep5] = float(col4[col3.index(nmid)])
    	dicip1["b"]["tail"][onsep1, onsep5] = float(col4[col3.index(ntail)])
    	print(col3.index(nhead),col3.index(nmid),col3.index(ntail),float(col4[col3.index(nhead)]),float(col4[col3.index(nmid)]),float(col4[col3.index(ntail)]))
    	f.close()

# Horizontal orbit offset for IP5
    else:
    	path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/hoff_f.IP5".format(pr = process,st = stg,num = n)
    	f = open(path, 'r')

    	col17, col18 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1])
    	col17 = list(col17)
    	col18 = list(col18)
    	dicip5["f"]["head"][onsep1, onsep5] = float(col18[col17.index(nhead)])
    	dicip5["f"]["medium"][onsep1, onsep5] = float(col18[col17.index(nmid)])
    	dicip5["f"]["tail"][onsep1, onsep5] = float(col18[col17.index(ntail)])
    	f.close()

    	path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/hoff_b.IP5".format(pr = process,st = stg,num=n)
    	f = open(path, 'r')

    	col19, col20 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1])
    	col19 = list(col19)
    	col20 = list(col20)
    	dicip5["b"]["head"][onsep1, onsep5] = float(col20[col19.index(nhead)])
    	dicip5["b"]["medium"][onsep1, onsep5] = float(col20[col19.index(nmid)])
    	dicip5["b"]["tail"][onsep1, onsep5] = float(col20[col19.index(ntail)])
    	f.close()


# Tunes and chromaticities
    path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/tune_f.list".format(pr = process,st = stg,num = n)
    f = open(path,'r')

    col5, col6, col7, col8, col9, col10 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1, 2, 3, 4, 5])
    col5 = list(col5)
    col6 = list(col6)
    col7 = list(col7)
    col8 = list(col8)
    col9 = list(col9)
    col10 = list(col10)

    tiph["f"]["head"][onsep1,onsep5] = float(col7[col5.index(nhead)])
    tiph["f"]["medium"][onsep1, onsep5] = float(col7[col5.index(nmid)])
    tiph["f"]["tail"][onsep1, onsep5] = float(col7[col5.index(ntail)])
    
    tipv["f"]["head"][onsep1,onsep5] = float(col8[col5.index(nhead)])
    tipv["f"]["medium"][onsep1, onsep5] = float(col8[col5.index(nmid)])
    tipv["f"]["tail"][onsep1, onsep5] = float(col8[col5.index(ntail)])
    
    
    ciph["f"]["head"][onsep1,onsep5] = float(col9[col5.index(nhead)])
    ciph["f"]["medium"][onsep1, onsep5] = float(col9[col5.index(nmid)])
    ciph["f"]["tail"][onsep1, onsep5] = float(col9[col5.index(ntail)])
    
    cipv["f"]["head"][onsep1,onsep5] = float(col10[col5.index(nhead)])
    cipv["f"]["medium"][onsep1, onsep5] = float(col10[col5.index(nmid)])
    cipv["f"]["tail"][onsep1, onsep5] = float(col10[col5.index(ntail)])
    

    f.close()
    
    path = "RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/tune_b.list".format(pr = process,st = stg,num = n)
    f = open(path,'r')

    col11, col12, col13, col14, col15, col16 = np.loadtxt(f, unpack=True, skiprows=0, usecols=[0, 1, 2, 3, 4, 5])
    col11 = list(col11)
    col12 = list(col12)
    col13 = list(col13)
    col14 = list(col14)
    col15 = list(col15)
    col16 = list(col16)

    tiph["b"]["head"][onsep1,onsep5] = float(col13[col11.index(nhead)])
    tiph["b"]["medium"][onsep1, onsep5] = float(col13[col11.index(nmid)])
    tiph["b"]["tail"][onsep1, onsep5] = float(col13[col11.index(ntail)])
    
    tipv["b"]["head"][onsep1,onsep5] = float(col14[col11.index(nhead)])
    tipv["b"]["medium"][onsep1, onsep5] = float(col14[col11.index(nmid)])
    tipv["b"]["tail"][onsep1, onsep5] = float(col14[col11.index(ntail)])
    
    
    ciph["b"]["head"][onsep1,onsep5] = float(col15[col11.index(nhead)])
    ciph["b"]["medium"][onsep1, onsep5] = float(col15[col11.index(nmid)])
    ciph["b"]["tail"][onsep1, onsep5] = float(col15[col11.index(ntail)])
    
    cipv["b"]["head"][onsep1,onsep5] = float(col16[col11.index(nhead)])
    cipv["b"]["medium"][onsep1, onsep5] = float(col16[col11.index(nmid)])
    cipv["b"]["tail"][onsep1, onsep5] = float(col16[col11.index(ntail)])
    

    f.close()




    subprocess.call(["rm","-r","RESULTS/MDcycle/onsep_{pr}_{st}_IP{num}/".format(pr = process,st = stg,num = n)])
    

    with open('dicip1.pickle', 'wb') as fileout:
        pk.dump(dicip1, fileout, protocol=pk.HIGHEST_PROTOCOL)
    
    with open('dicip5.pickle', 'wb') as fileout2:
    	pk.dump(dicip5, fileout2, protocol=pk.HIGHEST_PROTOCOL)

    with open('tiph.pickle', 'wb') as fileout3:
    	pk.dump(tiph, fileout3, protocol=pk.HIGHEST_PROTOCOL)
    
    with open('tipv.pickle', 'wb') as fileout4:
    	pk.dump(tipv, fileout4, protocol=pk.HIGHEST_PROTOCOL)
	
    with open('ciph.pickle', 'wb') as fileout5:
    	pk.dump(ciph, fileout5, protocol=pk.HIGHEST_PROTOCOL)
    
    with open('cipv.pickle', 'wb') as fileout6:
    	pk.dump(cipv, fileout6, protocol=pk.HIGHEST_PROTOCOL)
    

    output = '{st1} \t {st2} \t {st3} \n'.format(st1 = onsep1, st2 = onsep5, st3 = dicip1["f"]["head"][onsep1,onsep5])
    outfile.write(output)
    
def whole_matrix(onseprange,step,filename,dicip1,dicip5,prc,stg):

    for onsep1 in onseprange:
            subst = 'on_sep1              :=  {st} ;'.format(st = onsep1 - step)
            pattern = 'on_sep1              :=  {st} ;'.format(st = onsep1)

            replace(filename,subst,pattern)

            for onsep5 in onseprange:
                subst = 'on_sep5              :=  {st} ;'.format(st=onsep5 - step)
                pattern = 'on_sep5              :=  {st} ;'.format(st=onsep5)

                replace(filename, subst, pattern)
		
		f = open(filename,'r')
		content = f.read()
		print(content)
		f.close()

                dicip1,dicip5 = execute_train(dicip1,dicip5,onsep1,onsep5,prc,stg)


            subst5 = 'on_sep5              :=  {st} ;'.format(st=stop - step)
            pattern5 = 'on_sep5              :=  {st} ;'.format(st=start)

            replace(filename, subst5, pattern5)

    return dicip1,dicip5

def diagonal_matrix(n,outfile,onseprange,start,step,filename,dicip1,dicip5,tiph,tipv,ciph,cipv,prc,stg):
    
    for onsep in onseprange:
    	    subst = 'on_sep1              :=  {st} ;'.format(st = onsep - step)
            pattern = 'on_sep1              :=  {st} ;'.format(st = onsep)

            replace(filename,subst,pattern)
    	    
	    subst = 'on_sep5              :=  {st} ;'.format(st = onsep - step)
            pattern = 'on_sep5              :=  {st} ;'.format(st = onsep)

            replace(filename,subst,pattern)
            
            f = open(filename,'r')
	    content = f.read()
	    print('CONTENT')
	    print(content)
            f.close()

            onsepm = int(round((onsep - start - step )/step))

            execute_train(n,outfile,dicip1,dicip5,tiph,tipv,ciph,cipv,onsepm,onsepm,prc,stg)
    

    


def initialize_dics(stop):

    mtrx1 = np.zeros((stop, stop))
    mtrx2 = np.zeros((stop, stop))
    mtrx3 = np.zeros((stop, stop))
    mtrx4 = np.zeros((stop, stop))
    mtrx5 = np.zeros((stop, stop))
    mtrx6 = np.zeros((stop, stop))
    mtrx7 = np.zeros((stop, stop))
    mtrx8 = np.zeros((stop, stop))
    mtrx9 = np.zeros((stop, stop))
    mtrx10 = np.zeros((stop, stop))
    mtrx11 = np.zeros((stop, stop))
    mtrx12 = np.zeros((stop, stop))

    dicip1 = {"f": {"head": mtrx1, "tail": mtrx2, "medium": mtrx3}, "b": {"head": mtrx4, "tail": mtrx5, "medium": mtrx6}}
    dicip5 = {"f": {"head": mtrx7, "tail": mtrx8, "medium": mtrx9}, "b": {"head": mtrx10, "tail": mtrx11, "medium": mtrx12}}

    return dicip1,dicip5


def main():

    init = False
    prc = "collision"
    stg = "std"

    step = 0.1
    start = -6.0 - step
    stop = 6.0 + step
    n = 1
    

    filename = "MAD_PART/MDcycle/on_sep.{pr}".format(pr = prc)

    dicip1,dicip5 = initialize_dics(int((stop-start)/step))
    tiph,tipv = initialize_dics(int((stop-start)/step))
    ciph,cipv = initialize_dics(int((stop-start)/step))  

    if init == True:
        initialize(filename, start, -6.0)

    else:
        #onseprange = np.arange(start + step,stop,step)
	
	onsep = np.arange(start + step,-step,step)
	onsep2 = np.arange(-step,stop,step)


	onsep = list(onsep) + list(onsep2)
	onsep = [ '%.2f' % elem for elem in onsep ]
	onseprange = [ float(elem) for elem in onsep ]
	
	for i in onseprange:
		print(int(round((i - start - step )/step)))

        #dicip1,dicip5 = whole_matrix(matrixrange,filename,dicip1,dicip5,prc,stg)
	filename2 = 'follow.dat'
	outfile = open(filename2,'w')

	diagonal_matrix(n,outfile,onseprange,start,step,filename,dicip1,dicip5,tiph,tipv,ciph,cipv,prc,stg)
	outfile.close()

        

    
if __name__ == "__main__":
    main()
