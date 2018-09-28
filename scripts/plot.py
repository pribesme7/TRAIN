'''
Created on Apr 4, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import sys
import os

def main(arg1):
	print('Number of arguments:', len(sys.argv), 'arguments.')
	version = arg1
	resultFolder = arg1
	print("Plotting from the : "+resultFolder);

	if os.path.isfile(resultFolder + 'lumi.IP1'):
		slot1, lumiIP1 = np.loadtxt(resultFolder + 'lumi.IP1', unpack=True, skiprows=0, usecols=[0, 1])
	if os.path.isfile(resultFolder + 'lumi.IP5'):
		slot1, lumiIP5 = np.loadtxt(resultFolder + 'lumi.IP5', unpack=True, skiprows=0, usecols=[0, 1])
	 
	# ## Horizontal TUNE

	slot1, tune1f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 3])
	slot1, tune4f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 4])
	slot1, tune5f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 5])

	slot1, tune1b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 3])
	slot1, tune4b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 4])
	slot1, tune5b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 5])


	plt.figure(1);
	plt.subplot(221);
	plt.plot(slot1, tune2f, 'k.');
	plt.plot(slot1, tune3f, 'r.');
	plt.xlabel('slot number (25 ns)');
	plt.grid(True);
	plt.ylabel('  ');
	plt.legend(["tuneH B1" , "tuneV B1"])
	ax = plt.gca()

	plt.subplot(222);
	plt.plot(slot1, tune4f, 'k.');
	plt.plot(slot1, tune5f, 'r.');
	plt.grid(True);
	plt.xlabel('slot number (25 ns)');
	plt.ylabel('  ');
	plt.legend(["chromaX B1", "chromaV B1" ])
	ax = plt.gca()

	plt.subplot(223);
	plt.plot(slot1, tune2b, 'k.');
	plt.plot(slot1, tune3b, 'r.');
	plt.grid(True);
	plt.xlabel('slot number (25 ns)');
	plt.ylabel(' tune  ');
	plt.legend(["tuneH B2", "tuneV B2"])
	ax = plt.gca()

	plt.subplot(224);
	plt.plot(slot1, tune4b, 'k.');
	plt.plot(slot1, tune5b, 'r.');
	plt.xlabel('slot number (25 ns)');
	plt.grid(True);
	plt.ylabel('  ');
	plt.legend(["chromaX B2", "chromaV B2"])
	ax = plt.gca()


	plt.figure(5);
	plt.subplot(211);
	n, bins, patches = plt.hist(tune2f, 50, normed=0, facecolor='b', alpha=0.25);
	n, bins, patches = plt.hist(tune2b, 50, normed=0, facecolor='r', alpha=0.25);
	plt.ylabel('counts');
	plt.xlabel('Horizontal tune ');
	plt.legend(loc='upper right')
	plt.title(version);
	plt.grid(True);
	plt.subplot(212);
	plt.grid(True);
	n, bins, patches = plt.hist(tune3f, 50, normed=0, facecolor='b', alpha=0.75);
	n, bins, patches = plt.hist(tune3b, 50, normed=0, facecolor='r', alpha=0.75);
	plt.xlabel('Vertical tune ');
	plt.ylabel('counts');


	plt.figure(3);
	plt.subplot(211);
	plt.plot(slot1, tune1f, 'k.');
	plt.xlabel('slot number (25 ns)');
	plt.ylabel('  ');
	plt.legend(["bunch current"])
	ax = plt.gca()
	plt.subplot(212);
	if os.path.isfile(resultFolder + 'lumi.IP1'):
		plt.plot(slot1, lumiIP1, 'k.');
	if os.path.isfile(resultFolder + 'lumi.IP5'):
		plt.plot(slot1, lumiIP5, 'r.');
	plt.xlabel('slot number (25 ns)');
	plt.ylabel('% of lumi ');
	if (os.path.isfile(resultFolder + 'lumi.IP1') and os.path.isfile(resultFolder + 'lumi.IP5')) :
		plt.legend(["lumi IP1", "lumi IP5"])
	ax = plt.gca()
	plt.show();

if __name__ == "__main__":
    main(sys.argv[1])





