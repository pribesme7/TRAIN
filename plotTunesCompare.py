'''
Created on Apr 4, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import sys
import os

def main(resultFolder, resultFolder2):

	version = resultFolder +" vs " + resultFolder2

	print("Plotting from the : "+resultFolder + " and " + resultFolder2);
	 
	# ## Horizontal TUNE

	slot1, tune1f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3f = np.loadtxt(resultFolder + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 3])

	slot1, tune1b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3b = np.loadtxt(resultFolder + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 3])

	slot1, tune1f2 = np.loadtxt(resultFolder2 + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2f2 = np.loadtxt(resultFolder2 + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3f2 = np.loadtxt(resultFolder2 + 'tune_f.list', unpack=True, skiprows=0, usecols=[0, 3])

	slot1, tune1b2 = np.loadtxt(resultFolder2 + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 1])
	slot1, tune2b2 = np.loadtxt(resultFolder2 + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 2])
	slot1, tune3b2 = np.loadtxt(resultFolder2 + 'tune_b.list', unpack=True, skiprows=0, usecols=[0, 3])


	plt.figure(5);
	plt.subplot(211);
	n, bins, patches = plt.hist(tune2f, 50, normed=0, facecolor='b', alpha=0.75, label=" B1 from "+resultFolder);
	n, bins, patches = plt.hist(tune2b, 50, normed=0, facecolor='r', alpha=0.75, label=" B2 from "+resultFolder);
	n, bins, patches = plt.hist(tune2f2, 50, normed=0, facecolor='b', alpha=0.5, label=" B1 from "+resultFolder2);
	n, bins, patches = plt.hist(tune2b2, 50, normed=0, facecolor='r', alpha=0.5, label=" B2 from "+resultFolder2);
#	plt.plot((0.31, 0.31), (0, 2000), 'k-')
	plt.legend(loc='upper right')
	plt.xlabel('horizontal tune ');
	plt.ylabel('counts');
	plt.grid(True);

	plt.subplot(212);
	n, bins, patches = plt.hist(tune3f2, 50, normed=0, facecolor='b', alpha=0.5, label=" B1 from "+resultFolder2);
	n, bins, patches = plt.hist(tune3b2, 50, normed=0, facecolor='r', alpha=0.5, label=" B2 from "+resultFolder);
	n, bins, patches = plt.hist(tune3f, 50, normed=0, facecolor='b', alpha=0.75, label=" B1 from "+resultFolder2);
	n, bins, patches = plt.hist(tune3b, 50, normed=0, facecolor='r', alpha=0.75, label=" B2 from "+resultFolder);
#	plt.plot((0.32, 0.32), (0, 2000), 'k-')
	plt.xlabel('vertical tune ');
	plt.ylabel('counts');
	plt.grid(True);

	ax = plt.gca()
	plt.show();


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])





