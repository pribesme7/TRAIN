'''
Created on Apr 4, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import plotVerHorOffsets
import sys


def main(resultFolder, resultFolder2):
#	mainWithDisplay(resultFolder, resultFolder2, True, "dataFromOPScan4440");
	#same filling scheme
	#mainWithDisplay(resultFolder, resultFolder2, True, "dataFromOPScan5181");
#	mainWithDisplay(resultFolder, resultFolder2, True, "dataFromOPScan_5173");

	mainWithDisplay(resultFolder, resultFolder2, True, "dataFromMD_5137_185");
#	mainWithDisplay(resultFolder, resultFolder2, True, "dataFromMD_5137_105");

def mainWithDataFile(resultFolder, resultFolder2, dataFileToDisplay):
	mainWithDisplay(resultFolder, resultFolder2, True, dataFileToDisplay);

def printFirst(array, length):
	index = 1;
	for one in array:
		if index == length:
			break
		print(one)
		index+=1 

def mainWithDisplay(resultFolder, resultFolder2, display, dataFileToDisplay):
	version = resultFolder+' compared to '+resultFolder2;
#	referenceSlotNumber = 1+12+20+20; # 4440
	referenceSlotNumber = 1+20 #+ 48 +20; # MD 5137
	beamsize = 0.015;
	amplFactor = 2;

	###############
	#### TRAIN data

	slot1, IP1hoffIP1sig = np.loadtxt(resultFolder + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP1voffIP1sig = np.loadtxt(resultFolder + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);

	slot12, IP1hoffIP1sig2 = np.loadtxt(resultFolder2 + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot22, IP1voffIP1sig2 = np.loadtxt(resultFolder2 + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	slot12, IP5hoffIP5sig2 = np.loadtxt(resultFolder2 + 'hsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot22, IP5voffIP5sig2 = np.loadtxt(resultFolder2 + 'vsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);

	vertical = np.loadtxt('RESULTS_COMPARE_NOMINAL/'+dataFileToDisplay+'_V.csv', unpack=True, delimiter=',', skiprows=0, usecols=[0]);
	verticalErr = np.loadtxt('RESULTS_COMPARE_NOMINAL/'+dataFileToDisplay+'_V.csv', unpack=True, delimiter=',', skiprows=0, usecols=[1]);  
	horizontal = np.loadtxt('RESULTS_COMPARE_NOMINAL/'+dataFileToDisplay+'_H.csv', unpack=True, delimiter=',', skiprows=0, usecols=[0]);
	horizontalErr = np.loadtxt('RESULTS_COMPARE_NOMINAL/'+dataFileToDisplay+'_H.csv', unpack=True, delimiter=',', skiprows=0, usecols=[1]);  
		
	#fix from matlab export file (all slots) to only used slots
	validHorizontal = []; 
	validHorizontalErr= [];
	slotOptima = [];

	slotIndex = 0
	for one in horizontal:
		if one != 0.0:
			validHorizontal.append(amplFactor*(one/beamsize))
		slotOptima.append(slotIndex-1)
		slotIndex = slotIndex + 1			

	for one in horizontalErr:
		if one != 0.0:
			validHorizontalErr.append((one/beamsize))

	for x in range (1,13):
		validHorizontal.append(0)
		slotOptima.append(slotIndex-1)
		validHorizontalErr.append(0)
		slotIndex = slotIndex + 1			
	
	#fix from matlab export file (all slots) to only used slots
	validVertical = [];
	validVerticalErr = [];
	slotOptima = [];
	
	slotIndex = 0
	
	for one in vertical:
		if one != 0.0:
			slotOptima.append(slotIndex-1)
			validVertical.append(amplFactor * (one/beamsize))
		slotIndex = slotIndex + 1	

	for one in verticalErr:
		if one != 0.0:
			validVerticalErr.append(amplFactor * (one/beamsize))


	for x in range (1,13):
		slotOptima.append(slotIndex-1)
		validVertical.append(0)
		validVerticalErr.append(0)	
		slotIndex = slotIndex + 1	

	
	print(slot1)
	print(slot12)
	print(slotOptima)

	font = {'family' : 'serif',
        'size'   : 26}
######  SEPARATION PLANES
	plt.figure(2,figsize=(22,12));	
	plt.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.08);	
	plt.rc('font', **font);
	ax1 = plt.subplot(211);
#	plt.title('Horizontal separation at IP1 (Sep plane) \n' + version);
	plt.title('Horizontal separation at IP1 (Sep plane)');
	ax1.plot(slot1, np.multiply(plotVerHorOffsets.normalizeArray(IP1hoffIP1sig, referenceSlotNumber),1), 'k.');
	ax1.plot(slot12, np.multiply(plotVerHorOffsets.normalizeArray(IP1hoffIP1sig2, referenceSlotNumber),1), 'r.');
	#ax1.errorbar(slotOptima, plotVerHorOffsets.normalizeArray(validHorizontal, referenceSlotNumber), yerr=validHorizontalErr, color='g', fmt='o', ecolor='g');
	plt.ylabel('offset [$\sigma$]');
	plt.grid(True);
	plt.legend([resultFolder, resultFolder2]);
	ax = plt.gca()

	ax5 = plt.subplot(212, sharex=ax1);
#	plt.title('Vertical separation at IP5 (Sep plane) \n ' + version);
	plt.title('Vertical separation at IP5 (Sep plane)');
	ax5.plot(slot1, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig, referenceSlotNumber),  -1), 'k.');
	ax5.plot(slot12, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig2, referenceSlotNumber), -1), 'r.');
	ax5.errorbar(slotOptima, plotVerHorOffsets.normalizeArray(validVertical, referenceSlotNumber), yerr=validVerticalErr, color='g', fmt='o', ecolor='g');
	plt.grid(True);
	plt.xlabel('bunch position');
	plt.ylabel('offset [$\sigma$]');
	ax = plt.gca();
	plt.grid(True);
	
######  XING PLANES

	plt.figure(3,figsize=(22,12));
	plt.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.08);		
	plt.rc('font', **font);	
	ax2 = plt.subplot(211, sharex=ax1);
#	plt.title('Horizontal separation at IP5 (Xing plane) \n' + version);
	plt.title('Horizontal separation at IP5 (Xing plane)');
	ax2.plot(slot1, plotVerHorOffsets.normalizeArray(-IP5hoffIP5sig, referenceSlotNumber), 'k.');
	ax2.plot(slot12, plotVerHorOffsets.normalizeArray(-IP5hoffIP5sig2, referenceSlotNumber), 'r.');
	ax2.errorbar(slotOptima, np.multiply(plotVerHorOffsets.normalizeArray(validHorizontal, referenceSlotNumber),-1), yerr=validHorizontalErr, color='g', fmt='o', ecolor='g');
	plt.grid(True);
	plt.ylabel('offset [$\sigma$]');
	plt.legend([resultFolder, resultFolder2,"OP SCAN "+dataFileToDisplay]);
	ax = plt.gca();

	ax4 = plt.subplot(212, sharex=ax1);
#	plt.title('Vertical separation at IP1 (Xing plane) \n' + version);
	plt.title('Vertical separation at IP1 (Xing plane)');
	ax4.plot(slot1, plotVerHorOffsets.normalizeArray(-IP1voffIP1sig, referenceSlotNumber), 'k.');
	ax4.plot(slot12, plotVerHorOffsets.normalizeArray(-IP1voffIP1sig2, referenceSlotNumber), 'r.');
	plt.grid(True);
	plt.xlabel('bunch position [1]');
	plt.ylabel('offset [$\sigma$]');
#	plt.legend([resultFolder, resultFolder2]);
	ax = plt.gca()

#for papers plots!
######## SINGLE PLOT SEPARATION PLANE
	plt.figure(4,figsize=(22,5));	
	plt.subplots_adjust(left=0.08, right=0.97, top=0.90, bottom=0.20);	
	plt.rc('font', **font);
#	plt.title('Vertical separation at IP5 (Sep plane) \n ' + version);
	plt.title('Vertical separation at IP5 (Separation plane)');
	plt.plot(slot1, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig, referenceSlotNumber),  -1), 'k.');
	plt.plot(slot12, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig2, referenceSlotNumber), -1), 'r.');
	plt.errorbar(slotOptima, plotVerHorOffsets.normalizeArray(validVertical, referenceSlotNumber), yerr=validVerticalErr, color='g', fmt='o', ecolor='g');
	plt.grid(True);
	plt.xlabel('bunch position [1]');
	plt.ylabel('offset [$\sigma$]');
	plt.grid(True);
	plt.legend(["Simulation 185$\mu$rad", "Simulation 105$\mu$rad" ,"Emittance scan "+  dataFileToDisplay]);

######## IP5 PLOTS, SEP and XING.
	plt.figure(5,figsize=(22,12));
	plt.subplots_adjust(left=0.08, right=0.97, top=0.94, bottom=0.08);		
	plt.rc('font', **font);	
	ax66 = plt.subplot(211, sharex=ax1);
	plt.title('Horizontal separation at IP5 (Xing plane)');
	ax66.plot(slot1, plotVerHorOffsets.normalizeArray(-IP5hoffIP5sig, referenceSlotNumber), 'k.');
	ax66.plot(slot12, plotVerHorOffsets.normalizeArray(-IP5hoffIP5sig2, referenceSlotNumber), 'r.');
	ax66.errorbar(slotOptima, np.multiply(plotVerHorOffsets.normalizeArray(validHorizontal, referenceSlotNumber),-1), yerr=validHorizontalErr, color='g', fmt='o', ecolor='g');
	plt.grid(True);
	plt.ylabel('offset [$\sigma$]');
	plt.legend(["Simulation IP1 and IP5 incl.", "Simulation IP1, IP5 and IP8 incl." ,"Emittance Scan "+dataFileToDisplay]);
	ax77 = plt.subplot(212, sharex=ax1);
	plt.title('Vertical separation at IP5 (Separation plane)');
	ax77.plot(slot1, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig, referenceSlotNumber),  -1), 'k.');
	ax77.plot(slot12, np.multiply(plotVerHorOffsets.normalizeArray(IP5voffIP5sig2, referenceSlotNumber), -1), 'r.');
	ax77.errorbar(slotOptima, plotVerHorOffsets.normalizeArray(validVertical, referenceSlotNumber), yerr=validVerticalErr, color='g', fmt='o', ecolor='g');
	plt.grid(True);
	plt.xlabel('bunch position [1]');
	plt.ylabel('offset [$\sigma$]');
	ax = plt.gca();
	plt.grid(True);
	ax = plt.gca()
	if display:
		plt.show();


if __name__ == "__main__":
	if len(sys.argv) == 4:
		mainWithDisplay(sys.argv[1], sys.argv[2], False);
	else:
		main(sys.argv[1], sys.argv[2])
