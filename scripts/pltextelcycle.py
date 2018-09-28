import numpy as np
import matplotlib.pyplot as plt
import pickle as p
from matplotlib import colors as mcolors
from matplotlib.pyplot import xticks
from pylab import rcParams
rcParams['figure.figsize'] = 15, 5
axes = plt.gca()
axes.get_yaxis().get_major_formatter().set_scientific(True)
pdict = p.load(open( "save.p", "rb" ))
print(pdict)
colors = ['aqua','blue','coral','goldenrod','lavender','lime','maroon','red','magenta','teal','indigo','orchid','sienna']
elements = pdict['B1']['std'].keys()
fillsch = ['8','B','std']
resultFolder = 'Extraelements/'
ip = 'IP15'
print(elements)
for k in ['B1','B2']:
	for j in fillsch: 
		for l in ['h','v']:
 			fig  = plt.figure()
			for i in range(len(elements)):
				plt.plot(range(len(pdict[k][j][elements[i]][l])),pdict[k][j][elements[i]][l],'o',color = colors[i],label=elements[i],markersize=12)
				plt.plot(range(len(pdict[k][j][elements[i]][l])),pdict[k][j][elements[i]][l],color = colors[i],markersize=10)
			plt.legend()
			locs, labels = xticks()
			xticks(np.arange(5), ('Injection','Flat Top','Collision (64cm)','Squeeze','Physics'))
			plt.ylabel(r'Pick to pick orbit offset (um)',fontsize=20)
			plt.xticks(fontsize=18, rotation=0)
			plt.yticks(fontsize=18, rotation=0)
                        if l == 'h':
				strtitle = 'horizontal'
     			else:
				strtitle = 'vertical'
				
              		plt.title('{st1}, {st2} plane'.format(st1=k, st2=strtitle),fontsize=22)
			plt.grid(True)
#plt.legend(fontsize=13)
			fig.savefig(resultFolder + 'cycle_{st1}_{st2}_{st3}'.format(st1=k,st2=j,st3=l)+ip+ '.png')
			plt.show()


