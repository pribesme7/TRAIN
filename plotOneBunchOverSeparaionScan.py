'''
Created on Apr 12, 2017

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import plotVerHorOffsets
import sys, os

def plotFromResultsFolderForBunchId(resultsFolder, bunchId, ipToPlot = 'IP5'):


##### analytical formula
	xSep, coherent = np.loadtxt('data.out',delimiter=',', unpack=True, skiprows=0, usecols=[0, 7]);
	xSepIn, incoherent = np.loadtxt('data.out',delimiter=',', unpack=True, skiprows=0, usecols=[0, 2]);	


	beamSize = 0.000015;
	separation = [];
	separationSimulated = [];	
	b1OffsetV = [];
	b2OffsetV = [];
	b1OffsetH = [];
	b2OffsetH = [];
	
	for root, dirs, files in os.walk(resultsFolder):
		for directory in sorted(dirs):
			scanSeparation = int(directory.split(".")[-1]);
		
			print(directory)
			print (scanSeparation,'[um], beam size->', scanSeparation*1e-6/beamSize);
			separation.append(scanSeparation);

#			slot1, IP1hoffIP1sig = np.loadtxt(resultsFolder + directory + '/hsep_mu.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
#			slot2, IP1voffIP1sig = np.loadtxt(resultsFolder + directory + '/vsep_mu.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
			
			slotH, hSepSig = np.loadtxt(resultsFolder + directory + '/hsep_sig.' + ipToPlot , unpack=True, skiprows=0, usecols=[0, 1]);
			slotV, vSepSig = np.loadtxt(resultsFolder + directory + '/vsep_sig.' + ipToPlot, unpack=True, skiprows=0, usecols=[0, 1]);


#			slotb1h, hOffB1 = np.loadtxt(resultsFolder + directory +  '/hoff_f.' + ipToPlot, unpack=True, skiprows=0, usecols=[0, 1]);
			slotb1h, vOffB1 = np.loadtxt(resultsFolder + directory +  '/voff_f.' + ipToPlot, unpack=True, skiprows=0, usecols=[0, 1]);

#			slotb2h, hOffB2 = np.loadtxt(resultsFolder + directory + '/hoff_b.' + ipToPlot, unpack=True, skiprows=0, usecols=[0, 1]);
			slotb2v, vOffB2 = np.loadtxt(resultsFolder + directory + '/voff_b.' + ipToPlot, unpack=True, skiprows=0, usecols=[0, 1]);


			separationSimulated.append(hSepSig[bunchId])		
			b1OffsetV.append(vOffB1[bunchId])
			b2OffsetV.append(vOffB2[bunchId])


#	print(separation);
#	print(separationSimulated);	
#	print(b1OffsetV);	
#	print(b2OffsetV);	
#	print(b1OffsetH);	
#	print(b2OffsetH);	

	### pick one
#	separations = abs( np.array(separationSimulated) - separationSimulated[0] ); # reference it to the first step
	separations = np.array(b1OffsetV) - b1OffsetV[0] 			
#	separations = np.array(b2OffsetV) - b2OffsetV[0]
		
	## subtract the request separataion to see the rest		
	# factor 2 for total separation!!
	separationsNormalized = 2.0 * (separations - np.array(separation)) * 3.7;


############### PLOTING
			
	font = {'family' : 'serif',
        'size'   : 26}
	plt.figure(1,figsize=(15,6));	

	plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.15);	
	plt.rc('font', **font);	
	plt.title('Offsets in separation plane at %s'%ipToPlot);	
	plt.ylabel('BB offset [$\mu$m]');
	plt.xlabel('Separation [$\sigma$]');	#
	plt.plot(np.array(separation)/15.0, separationsNormalized, '-bx');
	plt.plot(np.array(xSep)/15.0, coherent, '-g');
#	plt.plot(np.array(xSepIn)/15.0, incoherent, '-r');	
	plt.xlim(0, 4);
	plt.ylim(0, 0.3);
	plt.legend(["TRAIN Coherent kick", "Analytical for Coherent kick", "Analytical for Incoherent kick"]);	
	plt.grid(True);
	
	
#	plt.figure(2,figsize=(15,6));	
#	plt.plot(separation, separations, 'r.');
#	plt.grid(True);
	
#	plt.figure(3,figsize=(15,6));	
#	plt.plot(xSep, coherent, '-g')
#	plt.grid(True);

	
#	plt.ylabel('offset [$\sigma$]');
#	plt.ylabel('offset [$\mu$m]');
#	plt.xlabel('separation [$\sigma$]');	
#	ax2 = plt.subplot(212, sharex=ax1);	
#	plt.show();




if __name__ == "__main__":
#		plotFromResultsFolderForBunchId(sys.argv[1], int(sys.argv[2]), ipToPlot='IP1')
		plotFromResultsFolderForBunchId(sys.argv[1], int(sys.argv[2]), ipToPlot='IP5')	
		
		plt.show();
