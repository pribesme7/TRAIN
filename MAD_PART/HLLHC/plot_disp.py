import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
rcParams['figure.figsize'] = 20, 5
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] = 16

path1 = 'collision64_train.optf.disp'
path2 = 'collision64_train.optb.disp'
path3 = 'collision64_train.optf.nodisp'
path4 = 'collision64_train.optb.nodisp'
path1 = 'train_stable_64_ondisp1.optf'
path2 = 'train_stable_64_ondisp1.optb'
path3 = 'train.optf'
path4 = 'train.optb'
sip1 = 1.999416239999800018e+04
sip2 = 2.332659898361884916e+04
sip5 = 6.664568432756300354e+03
sip8 = 1.665065818362084974e+04


s1d,x1d,dx1d,y1d,dy1d = np.loadtxt(path1,unpack=True,usecols=[1,2,5,6,9])
s2d,x2d,dx2d,y2d,dy2d = np.loadtxt(path2,unpack=True,usecols=[1,2,5,6,9])
s1,x1,dx1,y1,dy1 = np.loadtxt(path3,unpack=True,usecols=[1,2,5,6,9],skiprows=47)
s2,x2,dx2,y2,dy2 = np.loadtxt(path4,unpack=True,usecols=[1,2,5,6,9],skiprows=47)

plt.plot(s1d,x1d)
plt.plot(s1,x1)
plt.plot(s1d,y1d)
plt.plot(s1,y1)
plt.plot(s2d,x2d)
plt.plot(s2,x2)
plt.plot(s2d,y2d)
plt.plot(s2,y2)
plt.plot(sip1,0,'o',markersize=5)
plt.plot(sip2,0,'o',markersize=5)
plt.plot(sip5,0,'o',markersize=5)
plt.plot(sip8,0,'o',markersize=5)
plt.show()

plt.plot(s1d,dx1d,'--',label='on_disp = 1, B1')
plt.plot(s2d,dx2d,'--',label='on_disp = 1, B2')
plt.plot(s1,dx1,label='on_disp = 0, B1')
plt.plot(s2,dx2,label='on_disp = 0, B2')
plt.plot(sip1,0,'o',markersize=5)
plt.plot(sip2,0,'o',markersize=5)
plt.plot(sip5,0,'o',markersize=5)
plt.plot(sip8,0,'o',markersize=5)
plt.xlabel('s (m)',fontsize=18)
plt.grid(True)
plt.ylabel('Horizontal dispersion (m)',fontsize=18)
plt.legend()
plt.show()
plt.plot(s1d,dy1d,'--',label='on_disp = 1, B1')
plt.plot(s2d,dy2d,'--',label='on_disp = 1, B2')
plt.plot(s1,dy1,label='on_disp = 0, B1')
plt.plot(s2,dy2,label='on_disp = 0, B2')
plt.plot(sip1,0,'o',markersize=5)
plt.plot(sip2,0,'o',markersize=5)
plt.plot(sip5,0,'o',markersize=5)
plt.plot(sip8,0,'o',markersize=5)
plt.xlabel('s (m)',fontsize=18)
plt.grid(True)
plt.ylabel('Vertical dispersion (m)',fontsize=18)
plt.legend()
plt.show()
