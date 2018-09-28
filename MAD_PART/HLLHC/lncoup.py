import numpy as np
import matplotlib.pyplot as plt
from collections import deque
def comp_dqmin(npart,gamma,betxf0,betyf0, betxb0,betyb0,xf0,xb0,yf0,yb0):
    r0 = 1.53469857 * 10 ** (-18)
    en = 2.5 * 10 ** (-6)


    sigmaxf = np.sqrt(en * betxf0 / gamma)
    sigmayf = np.sqrt(en * betyf0 / gamma)
    sigmaf = np.sqrt(sigmaxf ** 2 + sigmayf ** 2)

    sigmaxb = np.sqrt(en * betxb0 / gamma)
    sigmayb = np.sqrt(en * betyb0 / gamma)
    sigmab = np.sqrt(sigmaxb ** 2 + sigmayb ** 2)
    sigma = np.sqrt(sigmaf ** 2 + sigmab ** 2)
    d = np.sqrt((xf0 - xb0) ** 2 + (yf0 - yb0) ** 2)/sigma
    alpha = abs(np.arctan((yf0 - yb0) / (xf0 - xb0)))
    print('d,alpha',d,alpha*180/np.pi)
    dqmin =  2*r0 * npart / (np.pi * en) * np.sin(alpha) * np.cos(alpha) / d ** 2
    dqtrial = r0 * npart / (np.pi * en) * np.sin(alpha) * np.cos(alpha) / d ** 2
    print(dqmin,sigma*1000*20)
    return dqmin,d,dqtrial

fnb = open('dqminonsep_f.dat','r')
bnb = open('dqminonsep_b.dat','r')
dnb,x0fnb,y0fnb,dqminfnb,dqantnb,betxfnb,betyfnb,gammafnb,rf11nb,rf12nb,rf21nb,rf22nb,qxfnb,qyfnb = np.loadtxt(fnb,unpack=True, skiprows=0, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13])
dnb,x0bnb,y0bnb,dqminbnb,dqantnb,betxbnb,betybnb,gammabnb,rb11nb,rb12nb,rb21nb,rb22nb,qxbnb,qybnb = np.loadtxt(bnb,unpack=True, skiprows=0, usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13])
fnb.close()
bnb.close()

npart =  2E14

newdqan = []
dn = []
detf = []
detb = []
dqtrial = []
for i in range(len(dnb)):
    dq,d,dqt = comp_dqmin(npart,gammafnb[i],0.64,0.64, 0.64,0.64,x0fnb[i],x0bnb[i],y0fnb[i],y0bnb[i])
    newdqan.append(dq)
    dn.append(d)
    dqtrial.append(dqt)
axes = plt.gca()
axes.get_yaxis().get_major_formatter().set_scientific(True)
print(dqminfnb)
dqminfnb = np.array(list([0]) + list(dqminfnb[int(round(len(dn)/2)):-1]))
dqminbnb = np.array(list([0]) + list(dqminbnb[int(round(len(dn)/2)):-1]))
#dqminfnb = dqminfnb[int(round(len(dn)/2)):]
#dqminbnb = dqminbnb[int(round(len(dn)/2)):]
print(dqminfnb)
print(dn[int(round(len(dn)/2)):])
plt.plot(dn[int(round(len(dn)/2)):],dqminfnb,'bo',label = 'B1')
plt.plot(dn[int(round(len(dn)/2)):],dqminbnb,'ro',label = 'B2')
plt.plot(dn[int(round(len(dn)/2)):],newdqan[int(round(len(dn)/2)):],'g+',label = 'Analytic')
#plt.plot(dn[int(round(len(dn)/2)):],dqtrial[int(round(len(dn)/2)):],'y*',label = 'Analytic')
#axes.set_xlim([400,1400])
#axes.set_ylim([0.0,0.0003])
plt.legend()
plt.show()
