import numpy as np
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
	referenceSlotNumber = -1;
	slot1 = []
	slot11 = []
	slot2 = []
	slot21 = []
	IP5hoffB1 = []
	IP5voffB1 = []
	IP5hoffB2 = []
	IP5voffB2 = []
	IP5hoffIP5sig = []
	IP5voffIP5sig = []
	if os.path.isfile(resultFolder + 'hoff_f.IP5'):
		slot1, IP5hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP5voffB1 = np.loadtxt(resultFolder + 'voff_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP5voffB2 = np.loadtxt(resultFolder + 'voff_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP5', unpack=True, skiprows=0, usecols=[0, 1]);

	elif os.path.isfile(resultFolder + 'hoff_f.IP1'):
		slot1, IP5hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP5voffB1 = np.loadtxt(resultFolder + 'voff_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP5voffB2 = np.loadtxt(resultFolder + 'voff_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

        elif os.path.isfile(resultFolder + 'hoff_f.IP2'):
		slot1, IP5hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP5voffB1 = np.loadtxt(resultFolder + 'voff_f.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP5voffB2 = np.loadtxt(resultFolder + 'voff_b.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP2', unpack=True, skiprows=0, usecols=[0, 1]);
		print(len(slot1),len(IP5hoffB1))

	elif os.path.isfile(resultFolder + 'hoff_f.IP8'):
		slot1, IP5hoffB1 = np.loadtxt(resultFolder + 'hoff_f.IP8', unpack=True, skiprows=0, usecols=[0, 1]);
		slot11, IP5voffB1 = np.loadtxt(resultFolder + 'voff_f.IP8', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5hoffB2 = np.loadtxt(resultFolder + 'hoff_b.IP8', unpack=True, skiprows=0, usecols=[0, 1]);
		slot21, IP5voffB2 = np.loadtxt(resultFolder + 'voff_b.IP8', unpack=True, skiprows=0, usecols=[0, 1]);
		slot1, IP5hoffIP5sig = np.loadtxt(resultFolder + 'hsep_sig.IP8', unpack=True, skiprows=0, usecols=[0, 1]);
		slot2, IP5voffIP5sig = np.loadtxt(resultFolder + 'vsep_sig.IP8', unpack=True, skiprows=0, usecols=[0, 1]);

	plt.figure(1);

	ax3 = plt.subplot(221);
	ax3.plot(slot1, normalizeArray(IP5voffB1, referenceSlotNumber), 'k.');
	ax3.plot(slot1, normalizeArray(IP5voffB2, referenceSlotNumber), 'r.');
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
	plt.plot(slot1, normalizeArray(IP5hoffB1, referenceSlotNumber), 'k.');
	plt.plot(slot1, normalizeArray(IP5hoffB2, referenceSlotNumber), 'r.');
	plt.xlabel('slot id');
	plt.grid(True);
	plt.ylabel('um');
	plt.legend(["B1", "B2"]);
	plt.xlim([0.0, 1000.0]);
	ax = plt.gca();

	plt.figure(3);
	axx1 = plt.subplot(212, sharex=ax3);
	plt.title('Separation at IP5 ' + version);
	axx1.plot(slot1, normalizeArray(IP5hoffIP5sig, referenceSlotNumber), 'k.');
	axx1.plot(slot1, normalizeArray(IP5voffIP5sig, referenceSlotNumber), 'r.');
	plt.xlabel('slot id');
	plt.ylabel('sig');
	plt.legend(["Horizontal", "Vertical"]);
	plt.grid(True);
	ax = plt.gca();
	print('IN MYPLOT')
	plt.show()

if __name__ == "__main__":
    main(sys.argv[1])


