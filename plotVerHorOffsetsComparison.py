'''
Created on Apr 4, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import plotVerHorOffsets
import sys


def readOneDirectory(resultFolder):
	slot1, IP1hoffIP1sig = np.loadtxt(resultFolder + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot1, IP1voffIP1sig = np.loadtxt(resultFolder + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot1, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	return slot1, IP1hoffIP1sig, IP1voffIP1sig, IP5hoffIP5sig, IP5voffIP5sig


def main(resultFolder, resultFolder2):
	mainWithDisplay(resultFolder, resultFolder2, True);

def mainWithDisplay(resultFolder, resultFolder2, display):
        mainForMultipleFolders(display, (resultFolder, resultFolder2) )

def mainForMultipleFolders(display, folders):

	#referenceSlotNumber = 1+12+12+4*48+20;
	referenceSlotNumber = 1+12+20;
	plotData = [];
	foldersLabels = [];
	for oneFolder in folders:
	    plotData.append(readOneDirectory(oneFolder));
#	    foldersLabels.append(oneFolder.split("/")[1])
	    foldersLabels.append('$\epsilon$=' + oneFolder.split("/")[1].split('_')[-1].split("emitt")[-1] + '$\mu$mrad')

	numberOfPlots = len(folders);
#	print(range(0, numberOfPlots))

	font = {'family' : 'serif',
        'size'   : 24}

	plt.figure(2,figsize=(22,12));
	plt.subplots_adjust(left=0.08, right=0.97, top=0.9, bottom=0.1);	

	plt.rc('font', **font);	
	ax1 = plt.subplot(211);
	plt.title('Horizontal separation at IP1 (Sep plane)');
	for i in range(0, numberOfPlots):
		ax1.plot(plotData[i][0], plotVerHorOffsets.normalizeArray(-plotData[i][1], referenceSlotNumber),".");
	plt.grid(True);
	plt.ylabel('offset [$\sigma$]');
	plt.legend(list(foldersLabels), loc='upper right')
		
	ax = plt.gca()
	ax5 = plt.subplot(212, sharex=ax1);
	plt.title('Vertical separation at IP5 (Sep plane)');
	for i in range(0, numberOfPlots):
		ax5.plot(plotData[i][0], plotVerHorOffsets.normalizeArray(-plotData[i][4], referenceSlotNumber), '.');
	plt.grid(True);

	plt.ylabel('offset [$\sigma$]');
	plt.xlabel('slot id [1]');	
		

#	plt.legend(list(folders),loc='upper right', bbox_to_anchor=(0.5, -0.1))	
#	plt.subplots_adjust(left=0.04, right=0.97, top=0.94, bottom=0.3)
	
#	plt.savefig(resultFolder2+'IP1separations.png')
#	plt.savefig(resultFolder2+'IP1separations.pdf')

	plt.figure(3,figsize=(22,12));	
	plt.subplots_adjust(left=0.08, right=0.97, top=0.9, bottom=0.1);	
	plt.rc('font', **font);	

	
	ax4 = plt.subplot(211, sharex=ax1);	
	plt.title('Vertical separation at IP1 (Xing plane)');
	for i in range(0, numberOfPlots):
		ax4.plot(plotData[i][0], plotVerHorOffsets.normalizeArray(-plotData[i][2], referenceSlotNumber), '.');
	plt.grid(True);
	plt.ylabel('offset [$\sigma$]');

	ax2 = plt.subplot(212, sharex=ax1);
	
	plt.title('Horizontal separation at IP5 (Xing plane)');
	for i in range(0, numberOfPlots):
		ax2.plot(plotData[i][0], plotVerHorOffsets.normalizeArray(-plotData[i][3], referenceSlotNumber), '.');
	plt.grid(True);
	plt.ylabel('offset [$\sigma$]');


	plt.legend(list(foldersLabels),loc='upper right')
#	plt.subplots_adjust(left=0.04, right=0.97, top=0.94, bottom=0.3)

#	plt.savefig(resultFolder2+'IP5separations.png')
#	plt.savefig(resultFolder2+'IP5separations.pdf')
	
	
	if display:
		plt.show();


if __name__ == "__main__":
	
	if len(sys.argv) == 4:
		mainWithDisplay(sys.argv[1], sys.argv[2], False);
	else:
		main(sys.argv[1], sys.argv[2])
