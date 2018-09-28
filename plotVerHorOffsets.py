'''
Created on Apr 4, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import sys
import os.path

def normalizeArray(dataInSlots, slotId):
	if slotId == -1:
		return dataInSlots;
	normalized = [];
	for i in range(len(dataInSlots)):
		normalized.append(dataInSlots[i]-dataInSlots[slotId])
	return normalized;

def main(arg1):
	version = arg1
	resultFolder = arg1
	referenceSlotNumber = 1+12+12+15;
	referenceSlotNumber = 2+17;
	referenceSlotNumber = -1;

	if os.path.isfile(resultFolder + 'hoff_f.IP5'):
		slot1, IP5hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP5voffB1 = np.loadtxt(resultFolder + 'voff_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP5voffB2 = np.loadtxt(resultFolder + 'voff_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);

	if os.path.isfile(resultFolder + 'hoff_f.IP1'):
		slot1, IP1hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP1voffB1 = np.loadtxt(resultFolder + 'voff_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP1hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP1voffB2 = np.loadtxt(resultFolder + 'voff_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP1hoffIP1sig = np.loadtxt(resultFolder + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP1voffIP1sig = np.loadtxt(resultFolder + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	plt.figure(1);
	# IP5 is always in the results!
        ax3 = plt.subplot(221);
	if os.path.isfile(resultFolder + 'hoff_f.IP5'):
#		ax3 = plt.subplot(221);
		ax3.plot(slot11, normalizeArray(IP5voffB1, referenceSlotNumber), 'k.');
		ax3.plot(slot21, normalizeArray(IP5voffB2, referenceSlotNumber), 'r.');
		plt.title('Vertical offset at IP5 \n' + version);
		plt.xlabel('slot id');
		plt.ylabel('um');
		plt.grid(True);
		plt.legend(["B1", "B2"]);
		plt.xlim([0.0, 1000.0]);
		ax = plt.gca();
		####
		ax4 = plt.subplot(222, sharex=ax3);
		plt.title('Horizontal offset at IP5 \n' + version);
		plt.plot(slot11, normalizeArray(IP5hoffB1, referenceSlotNumber), 'k.');
		plt.plot(slot21, normalizeArray(IP5hoffB2, referenceSlotNumber), 'r.');
		plt.xlabel('slot id');
		plt.grid(True);
		plt.ylabel('um');
		plt.legend(["B1", "B2"]);
		plt.xlim([0.0, 1000.0]);
		ax = plt.gca();

	if os.path.isfile(resultFolder + 'hoff_f.IP1'):
		ax1 = plt.subplot(223, sharex=ax3);
		ax1.plot(slot11, normalizeArray(IP1voffB1, referenceSlotNumber), 'k.');
		ax1.plot(slot21, normalizeArray(IP1voffB2, referenceSlotNumber), 'r.');
		plt.title('Vertical offset at IP1 \n' + version)
		plt.xlabel('slot id');
		plt.ylabel('um');
		plt.grid(True);
		plt.legend(["B1", "B2"])
		ax = plt.gca()
		####
		ax2=plt.subplot(224, sharex=ax3);
		plt.title('Horizontal offset at IP1 \n' + version);
		ax2.plot(slot11, normalizeArray(IP1hoffB1, referenceSlotNumber), 'k.');
		ax2.plot(slot21, normalizeArray(IP1hoffB2, referenceSlotNumber), 'r.');
		plt.xlabel('slot id');
		plt.ylabel('um');
		plt.grid(True);
		plt.legend(["B1", "B2"]);
		#plt.xlim([0.0, 1000.0]);
		ax = plt.gca();

	plt.figure(3);
	# IP5 is always in the results!	
	if os.path.isfile(resultFolder + 'hoff_f.IP5'):
		axx1 = plt.subplot(212, sharex=ax3);
		plt.title('Separation at IP5 ' + version);
		axx1.plot(slot1, normalizeArray(-IP5hoffIP5sig, referenceSlotNumber), 'k.');
		axx1.plot(slot1, normalizeArray(-IP5voffIP5sig, referenceSlotNumber), 'r.');
		plt.xlabel('slot id');
		plt.ylabel('sig');
		plt.legend(["Horizontal", "Vertical"]);
		plt.grid(True);
		ax = plt.gca();	

	if os.path.isfile(resultFolder + 'hoff_f.IP1'):
		axx2 = plt.subplot(211, sharex=ax3);
		plt.title('Separation at IP1 ' + version);
		axx2.plot(slot1, normalizeArray(-IP1hoffIP1sig, referenceSlotNumber), 'k.');
		axx2.plot(slot1, normalizeArray(-IP1voffIP1sig, referenceSlotNumber), 'r.');
		plt.xlabel('slot id');
		plt.grid(True);
		plt.ylabel('sig');
		plt.legend(["Horizontal", "Vertical"]);
		ax = plt.gca()

	plt.show();


if __name__ == "__main__":
    main(sys.argv[1])


