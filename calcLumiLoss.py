'''
Created on Jun 21, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import plotVerHorOffsets
import sys
from math import exp
from math import sqrt

def beamSize(en):
	return sqrt(0.4 * en * 10e-6 / 6930) # beta 40cm, 6.5TeV

def normalizeArray(dataInSlots, slotId):
	if slotId == -1:
		return dataInSlots;
	normalized = [];
	for i in range(len(dataInSlots)):
		normalized.append(dataInSlots[i]-dataInSlots[slotId])
	return normalized;

def geometricFactor(separationX, separationY, crossingangle):
	return exp(- (separationX * separationX + separationY*separationY) / 4)

def getLumiPercentage(separationsSep, separationXing, slotId,  crossingangle):
	luminosityBB  = []
	nSeparationsSep=normalizeArray(separationsSep, slotId)
	nSeparationsXing=normalizeArray(separationXing, slotId)
	for bunch in range(len(nSeparationsSep)):
		luminosityBB.append(geometricFactor(nSeparationsSep[bunch], nSeparationsXing[bunch], crossingangle))
	return luminosityBB;

#def normalizeAndGetSeparation(one, two, slot):
#	beam1 = normalizeArray(one, slot)
#	beam2 = normalizeArray(two, slot)
#	offsets  = list(np.array(beam1)-np.array(beam2))
#	nSeparation = [];
#	for one in offsets:
#		nSeparation.append(abs(one)/beamSize(3.0))
#	return nSeparation;

def main(resultFolder, resultFolder2):
	
	crossingAngle = 185
	slotId  = 1 + 12 + 12 + 26; #somewhere in a middle of a first real train

# TRAIN sigmas of separation !! (??) what beam size!
	slot1, IP1hoffIP1sig = np.loadtxt(resultFolder + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP1voffIP1sig = np.loadtxt(resultFolder + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	slot1, IP1hoffIP1sig2 = np.loadtxt(resultFolder2 + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP1voffIP1sig2 = np.loadtxt(resultFolder2 + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

# TRAIN lumi loss
	slot1, lumiIP1 = np.loadtxt(resultFolder + 'lumi.IP5', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, lumiIP1f2 = np.loadtxt(resultFolder2 + 'lumi.IP5', unpack=True, skiprows=0, usecols=[0, 1])

# shifts of the centroid in um
#	slot1, IP1_H_f_disp = np.loadtxt(resultFolder + 'hdisp_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
#	slot2, IP1_H_b_disp = np.loadtxt(resultFolder + 'hdisp_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
#	slot1, IP1_V_f_disp = np.loadtxt(resultFolder + 'vdisp_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
#	slot2, IP1_V_b_disp = np.loadtxt(resultFolder + 'vdisp_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
# 	hNormalizedIP1 = normalizeAndGetSeparation(IP1_H_f_disp, IP1_H_b_disp, slotId)
# 	vNormalizedIP1 = normalizeAndGetSeparation(IP1_V_f_disp, IP1_V_b_disp, slotId)
#	lumiLossFromOffset = getLumiPercentage(hNormalizedIP1,  vNormalizedIP1  , slotId, crossingAngle)

	lumi  = getLumiPercentage(IP1hoffIP1sig,  IP1voffIP1sig  , slotId, crossingAngle)
	lumi2 = getLumiPercentage(IP1hoffIP1sig2, IP1voffIP1sig2 , slotId, crossingAngle)
	
	print( "Loss calculated from separations TRAIN")
	print("Average (for %d bunches) lumi loss %f. for dir1" % (len(lumi), sum(lumi)/float(len(lumi))))
	print("Average (for %d bunches) lumi loss %f. for dir2" % (len(lumi), sum(lumi2)/float(len(lumi2))))
	print("Loss from Lumi in TRAIN")
	print("Average (for %d bunches) lumi loss %f. for dir1" % (len(lumi), sum(lumiIP1)/float(len(lumiIP1))))
	print( "Average (for %d bunches) lumi loss %f. for dir2" % (len(lumi), sum(lumiIP1f2)/float(len(lumiIP1f2))))


#	plt.figure(1);
#	plt.title('Reduction per bunch slot');
#	plt.plot(slot1, lumi, 'b.');
#	#ax1.plot(slot1, lumiIP1, 'r.');
#	plt.plot(slot1, lumi2, 'r.');
#	#ax1.plot(slot1, lumiIP1f2, 'g.');
#	plt.legend([resultFolder,resultFolder2],loc='upper center', bbox_to_anchor=(0.5, -0.05))
##	plt.legend();
#	plt.xlabel('slot id');
#	plt.ylabel('%');
#	plt.grid(True);


#	plt.figure(2);
#	ax1 = plt.subplot(211);	
#	plt.title('Reduction per bunch slot');
#	ax1.plot(slot1, lumi, 'b.');
#	#ax1.plot(slot1, lumiIP1, 'r.');
#	ax1.plot(slot1, lumi2, 'r.');
#	#ax1.plot(slot1, lumiIP1f2, 'g.');
#	plt.legend([resultFolder,resultFolder2]);
#	plt.xlabel('slot id');
#	plt.ylabel('%');
#	plt.grid(True);
#	ax2 = plt.subplot(212, sharex=ax1);
#	ax2.plot(slot1, normalizeArray(IP1hoffIP1sig , slotId),'k.')
#	ax2.plot(slot1, normalizeArray(IP1voffIP1sig , slotId),'b.')
#	ax2.plot(slot1, normalizeArray(IP1hoffIP1sig2 , slotId),'g.')
#	ax2.plot(slot1, normalizeArray(IP1voffIP1sig2 , slotId),'r.')
#	plt.legend(["Dir1 H offset","Dir1 V offset","Dir2 H offset","Dir2 V offset"]);
	
#	plt.figure(3);
#	n, bins, patches = plt.hist(lumi, 50, normed=0, facecolor='b', alpha=0.75);
#	n, bins, patches = plt.hist(lumi2, 50, normed=0, facecolor='r', alpha=0.75);
#	plt.ylabel('bunches');
#	plt.xlabel('Lumi loss % ');
#	plt.legend(loc='upper right')
#	plt.grid(True);


	plt.figure(4,figsize=(18, 14));
	plt.title('Reduction histogram per bunch');
	ax1 = plt.subplot(211);	
	n, bins, patches = ax1.hist(lumi, 50, normed=0, facecolor='b', alpha=0.75);
	n, bins, patches = ax1.hist(lumi2, 50, normed=0, facecolor='r', alpha=0.75);
	plt.ylabel('bunches');
	plt.xlabel('Lumi loss % ');

	ax2 = plt.subplot(212);	
	plt.title('Reduction per bunch slot');
	ax2.plot(slot1, lumi, 'b.');
	ax2.plot(slot1, lumi2, 'r.');
	plt.xlabel('slot id');
	plt.ylabel('%');
	plt.grid(True);
	plt.legend([resultFolder,resultFolder2],loc='upper center', bbox_to_anchor=(0.5, -0.1))
	plt.subplots_adjust(left=0.1, right=0.95, top=0.96, bottom=0.14)
	plt.grid(True);	

	plt.show();
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])
