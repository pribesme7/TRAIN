'''
Created on Apr 24, 2016

@author: agorzaws
'''
#import tkFileDialog
from tkinter import filedialog
from tkinter import *
from multiprocessing import Process
import os
import plot
import plotVerHorOffsets
import plotVerHorOffsetsComparison
import plotVerHorOffsetsComparisonWithData
import plotTunesCompare
import plotVerHorControidPositionComparison
import calcLumiLoss

class OpenDialog(Frame):
    
	def __init__(self, master=None):
		self.root=master;
		self.createWidgets();

	def createWidgets(self):
		self.fFrame = Frame(self.root);
		self.bSelDir = Button(self.fFrame, text="load directory");
		self.bSelDir["command"] =  self.getDir;
		self.bSelDir.pack();

		self.bRefresh = Button(self.fFrame, text="refresh current directory");
		self.bRefresh["command"] =  self.refreshDir;
		self.bRefresh.pack();

		self.svDir = StringVar();
		global eDir;
		eDir  = Entry(self.fFrame, width=170, textvariable=self.svDir);
		eDir.pack();
        
		scrollbar = Scrollbar();
		scrollbar.pack(side=RIGHT, fill=Y);
		global lb1;
		lb1 = Listbox(self.fFrame,selectmode=MULTIPLE,yscrollcommand=scrollbar.set);
		lb1.pack(fill=BOTH, expand=1);
		lb1.config(yscrollcommand=scrollbar.set);
		lb1.pack();
		scrollbar.config(command=lb1.yview);

		self.bSelDir.pack();
		bTune = Button(self.fFrame, text="Load Tune view");
		bTune["command"] =  self.getTuneView;
		bTune.pack();

		bOrbit = Button(self.fFrame, text="Load Orbit view");
		bOrbit["command"] =  self.getOrbitView;
		bOrbit.pack();

		bLumi = Button(self.fFrame, text="Load orbit Luminosity view");
		bLumi["command"] =  self.getLumiView;
		bLumi.pack();

		bOrbitData = Button(self.fFrame, text="Load orbit comparison vs data");
		bOrbitData["command"] =  self.getOrbitData;
		bOrbitData.pack();

		self.fFrame.pack();

	def getDir(self):
		dir = filedialog.askdirectory();
		self.svDir.set(dir);
		self.loadFromDir(dir);

	def refreshDir(self):
		global eDir;
		self.loadFromDir(eDir.get());

	def loadFromDir(self, dirName):
		localDir = dirName.split('/')[-1];
		global lb1;
		lb1.delete(0, END);
		for root, dirs, files in os.walk(dirName):
			for dir in sorted(dirs):
				lb1.insert(END , localDir+"/"+dir+"/");
		lb1.pack();

	def getTuneView(self):
		global lb1;
		sel = lb1.curselection();
		if len(sel) == 1:
			actualDirToShow = lb1.get(sel);
			plot.main(actualDirToShow);
		elif len(sel) == 2:
			actualDirToShow = lb1.get(sel[0]);
			actualDirToShow2 = lb1.get(sel[1]);
			plotTunesCompare.main(actualDirToShow, actualDirToShow2);

	def getLumiView(self):
		global lb1;
		sel = lb1.curselection();
		if len(sel) == 2:
			actualDirToShow = lb1.get(sel[0]);		
			actualDirToShow2 = lb1.get(sel[1]);		
			calcLumiLoss.main(actualDirToShow, actualDirToShow2);

	def getOrbitData(self):
		global lb1;
		sel = lb1.curselection();

		if len(sel) == 2:
			actualDirToShow = lb1.get(sel[0]);
			actualDirToShow2 = lb1.get(sel[1]);		
			plotVerHorOffsetsComparisonWithData.main(actualDirToShow, actualDirToShow2);

	def getOrbitView(self):
		global lb1;
		sel = lb1.curselection();
		if len(sel) == 1:
			actualDirToShow = lb1.get(sel);
			p1 = Process(target = plotVerHorOffsets.main(actualDirToShow));
			p1.start();
			p2 = Process(target = plotVerHorControidPositionComparison.mainWithDisplay(actualDirToShow,"",True));
			p2.start();
	#		plotVerHorOffsets.main(actualDirToShow)
	#		plotVerHorControidPositionComparison.mainWithDisplay(actualDirToShow,"",True)
		else:
			print(sel);
			folders = [];
			for one in sel:
				folders.append(lb1.get(one));
			print(folders);
			plotVerHorOffsetsComparison.mainForMultipleFolders(True, folders);
	#	elif len(sel) == 2:
	#		actualDirToShow = lb1.get(sel[0])		
	#		actualDirToShow2 = lb1.get(sel[1])		
	#		plotVerHorOffsetsComparison.main(actualDirToShow, actualDirToShow2)

if __name__ == "__main__":
    root=Tk();
    od = OpenDialog(root);
    root.mainloop();
