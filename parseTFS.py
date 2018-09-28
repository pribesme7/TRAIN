import csv
import numpy as np
import matplotlib.pyplot as plt
#from plotLR import plotLR

def parseTFS(fileName):
    reader = csv.reader(open(fileName, 'r'), delimiter=' ');
    headers = 0;
    data = {};
    names = [];
    for row in reader:
        if len(row)>0:
            if len(row[0])>0:
                if row[0][0] == '@':
                    headers = headers + 1;
                elif row[0][0] == '*':
                    for i in np.arange(1,len(row)):
                        tmpName = row[i].strip();
                        if len(tmpName)>0:
                            names.append(tmpName);
                            data[tmpName] = [];
            else:
                dataCount = 0;
                for i in np.arange(len(row)):
                    tmpValue = row[i].strip();
                    if len(tmpValue) > 0:
                        try:
                            value = float(tmpValue);
                            data[names[dataCount]].append(value);
                        except ValueError:
                            data[names[dataCount]].append(tmpValue);
                        dataCount = dataCount + 1;
    return data;
             
        
if __name__ == "__main__":
    dataB1 = parseTFS('/home/xbuffat/MAD-X/thinTest/B1.twiss');
    dataB2 = parse('/home/xbuffat/MAD-X/thinTest/B2.twiss');
#    dataB1 = parse('/home/xbuffat/MAD-X/IP2/witnessBeam.tfs');
#    dataB2 = parse('/home/xbuffat/MAD-X/IP2/sourceBeam.tfs');
#    dataB1 = parse('/home/xbuffat/MAD-X/IP2/PHYSICS-TILTED-2012-XBtest_V1--B1--220.0--0.tfs');
#    dataB2 = parse('/home/xbuffat/MAD-X/IP2/PHYSICS-TILTED-2012-XBtest_V1--B2--220.0--0.tfs');
    
    print('parsing ok');
    
    nBPM = 0;
    for i in np.arange(len(dataB2['NAME'])):
        if 'BPM' in dataB2['NAME'][i]:
 #           print nBPM,dataB2['NAME'][i]
            nBPM = nBPM + 1;
    fig = plt.figure(0);
#    fig.add_subplot(211);
    plt.plot(dataB1['S'],dataB1['X'],label='x');
    plt.plot(dataB1['S'],dataB1['Y'],label='y');
    plt.ylabel('orbit B1 [m]');plt.legend(loc=0);
    #plt.twinx();
    #plt.plot(dataB1['S'],[180.0*np.arctan(dataB1['X'][k]/dataB1['Y'][k])/np.pi if np.abs(dataB1['Y'][k])>1E-6 else 0.0 for k in range(len(dataB1['X']))],label='angle');
    plt.xlabel('s [m]');plt.ylabel('xy angle');
    
    plt.figure(10);
    plt.plot(dataB1['X'],dataB1['Y']);
    plt.xlabel('x');plt.ylabel('y');
    
#    fig.add_subplot(212);
#    plt.plot(dataB2['S'],dataB2['X'],label='x');
#    plt.plot(dataB2['S'],dataB2['Y'],label='y');
#    plt.xlabel('s [m]');plt.ylabel('orbit B2 [m]');plt.legend(loc=0);
    fig = plt.figure(1);
#    fig.add_subplot(211);
    plt.plot(dataB1['S'],dataB1['BETX'],label='beta x');
    plt.plot(dataB1['S'],dataB1['BETY'],label='beta y');
    plt.xlabel('s [m]');plt.ylabel('beta B1 [m]');plt.legend(loc=0);
#    fig.add_subplot(212);
#    plt.plot(dataB2['S'],dataB2['BETX'],label='beta x');
#    plt.plot(dataB2['S'],dataB2['BETY'],label='beta y');
#    plt.xlabel('s [m]');plt.ylabel('beta B2 [m]');plt.legend(loc=0);

    gamma = 4000.0/0.938272013;
    emit = 2.5E-6/gamma;
    plt.figure(2);
    plt.plot(dataB1['S'],[np.sqrt((((dataB1['X'][k]-np.interp(dataB1['S'][k],dataB2['S'],dataB2['X']))/np.sqrt(np.interp(dataB1['S'][k],dataB2['S'],dataB2['BETX'])*emit))**2+((dataB1['Y'][k]-np.interp(dataB1['S'][k],dataB2['S'],dataB2['Y']))/np.sqrt(np.interp(dataB1['S'][k],dataB2['S'],dataB2['BETY'])*emit))**2)) for k in np.arange(len(dataB1['S']))],label='B1');
    plt.plot(dataB1['S'],[np.sqrt(((dataB1['X'][k]-np.interp(dataB1['S'][k],dataB2['S'],dataB2['X']))/np.sqrt(dataB1['BETX'][k]*emit))**2+((dataB1['Y'][k]-np.interp(dataB1['S'][k],dataB2['S'],dataB2['Y']))/np.sqrt(dataB1['BETY'][k]*emit))**2) for k in np.arange(len(dataB1['S']))],label='B2');
    plt.xlabel('S [m]');plt.ylabel('Separation [sigma]');plt.legend();
#    plotLR([0,30]);

    sIP8 = 23315.37898;
    #sIP8 = 3332.436584;
    dist = 1.0;
    XingHB1 = np.arctan((np.interp(sIP8+dist,dataB1['S'],dataB1['X'])-np.interp(sIP8-dist,dataB1['S'],dataB1['X']))/(2.0*dist))*1E6;
    XingVB1 = np.arctan((np.interp(sIP8+dist,dataB1['S'],dataB1['Y'])-np.interp(sIP8-dist,dataB1['S'],dataB1['Y']))/(2.0*dist))*1E6;
    XingHB2 = np.arctan((np.interp(sIP8+dist,dataB2['S'],dataB2['X'])-np.interp(sIP8-dist,dataB2['S'],dataB2['X']))/(2.0*dist))*1E6;
    XingVB2 = np.arctan((np.interp(sIP8+dist,dataB2['S'],dataB2['Y'])-np.interp(sIP8-dist,dataB2['S'],dataB2['Y']))/(2.0*dist))*1E6;
    print(XingHB1,XingVB1);
    print(XingHB2,XingVB2);
    
    plt.show();
