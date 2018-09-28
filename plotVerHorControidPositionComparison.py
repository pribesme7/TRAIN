'''
Created on Jun 14, 2016

@author: agorzaws
'''
import numpy as np
from numpy import sqrt, pi, exp, sign
import matplotlib.pyplot as plt
import plotVerHorOffsets
import sys


def main(resultFolder, resultFolder2):
	mainWithDisplay(resultFolder, resultFolder2, True);

def inverse(list):
	toReturn = [];
	for element in list:
		toReturn.append(-1*element)
	return toReturn;

def normalizeAndGetCentroid(one, two, slot):
	beam1 = plotVerHorOffsets.normalizeArray(one, slot)
	beam2 = plotVerHorOffsets.normalizeArray(two, slot)
	return list(np.array(beam1)-np.array(beam2))


def mainWithDisplay(resultFolder, resultFolder2, display):
	version = resultFolder;
	referenceSlotNumber = 1+12+12+15;
#	referenceSlotNumber = ;

	slot1, IP1_H_f_disp = np.loadtxt(resultFolder + 'hdisp_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP1_H_b_disp = np.loadtxt(resultFolder + 'hdisp_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot1, IP1_V_f_disp = np.loadtxt(resultFolder + 'vdisp_f.IP1', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP1_V_b_disp = np.loadtxt(resultFolder + 'vdisp_b.IP1', unpack=True, skiprows=0, usecols=[0, 1]);

	slot1, IP5_H_f_disp = np.loadtxt(resultFolder + 'hdisp_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP5_H_b_disp = np.loadtxt(resultFolder + 'hdisp_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot1, IP5_V_f_disp = np.loadtxt(resultFolder + 'vdisp_f.IP5', unpack=True, skiprows=0, usecols=[0, 1]);
	slot2, IP5_V_b_disp = np.loadtxt(resultFolder + 'vdisp_b.IP5', unpack=True, skiprows=0, usecols=[0, 1]);

	IP1_H_shift = normalizeAndGetCentroid(IP1_H_f_disp, IP1_H_b_disp, -1)
	IP1_V_shift = normalizeAndGetCentroid(IP1_V_f_disp, IP1_V_b_disp, referenceSlotNumber)

 	#IP1_H_shift = normalizeAndGetCentroid(IP1_H_f_disp, IP1_H_b_disp, -1)
 	#IP1_V_shift = normalizeAndGetCentroid(IP1_V_f_disp, IP1_V_b_disp, -1)
 	
	IP5_H_shift = normalizeAndGetCentroid(IP5_H_f_disp, IP5_H_b_disp, -1)
	IP5_V_shift = normalizeAndGetCentroid(IP5_V_f_disp, IP5_V_b_disp, referenceSlotNumber)

	plt.figure(1,figsize=(22,12));
	ax1 = plt.subplot(211);
	plt.title('Xing plane collision point shift \n' + version);
	ax1.plot(slot1, inverse(IP1_V_shift), 'k.',markersize=15);
	#ax1.plot(slot1, inverse(IP5_H_shift), 'r.',markersize=15);
	#ax1.legend(["IP1 (Vertical)", "IP5 (Horizontal)"]);
	plt.grid(True);

	ax2 = plt.subplot(212, sharex=ax1);
	plt.title('Separation plane collision point shift \n' + version);
	ax2.plot(slot1, inverse(IP1_H_shift), 'k.',markersize=15);
	#ax2.plot(slot1, inverse(IP5_V_shift), 'r.',markersize=15);
	#ax2.legend(["IP1 (Horizontal)", "IP5 (Vertical)"]);
	plt.grid(True);
	plt.ylabel('um');
	ax1 = plt.gca()


	if display:
		plt.show();


if __name__ == "__main__":
	
	if len(sys.argv) == 4:
		mainWithDisplay(sys.argv[1], sys.argv[2], False);
	else:
		main(sys.argv[1], sys.argv[2])
