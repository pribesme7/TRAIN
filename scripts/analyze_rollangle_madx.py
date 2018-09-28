import numpy as np

angle = np.loadtxt('MAD_PART/MDcycle/on_phi.collision41',unpack=True,usecols=[-2],skiprows=1)
s1,x1,y1 = np.loadtxt('MAD_PART/MDcycle/train.optf',unpack=True,usecols=[1,2,6],skiprows=47)
s2,x2,y2 = np.loadtxt('MAD_PART/MDcycle/train.optb',unpack=True,usecols=[1,2,6],skiprows=47)
s1 = list(s1)
s2 = list(s2)
sip1 = 1.999416239999800018e+04
sip5b2 = 6.664873167243699754e+03
sip5b1 = 6.664568432756300354e+03
sigma = np.sqrt(0.64*2.6*10**(-6)/7.460522527611191435e+03)
d1 = np.sqrt((x1[s1.index(sip1)-1]-x2[s2.index(sip1)-1])**2+(y1[s1.index(sip1)-1]-y2[s2.index(sip1)-1])**2)/sigma
alpha1 = np.arctan((y1[s1.index(sip1)-1]-y2[s2.index(sip1)-1])/x1[s1.index(sip1)-1]-x2[s2.index(sip1)-1])*180./np.pi

d5 = np.sqrt((x1[s1.index(sip5b1)-1]-x2[s2.index(sip5b2)-1])**2+(y1[s1.index(sip5b1)-1]-y2[s2.index(sip5b2)-1])**2)/sigma
alpha5 = np.arctan((y1[s1.index(sip5b1)-1]-y2[s2.index(sip5b2)-1])/x1[s1.index(sip5b1)-1]-x2[s2.index(sip5b2)-1])*180./np.pi

print(angle[0],angle[2],d1,alpha1,d5,alpha5)

d1 = np.sqrt((x1[s1.index(sip1)]-x2[s2.index(sip1)])**2+(y1[s1.index(sip1)]-y2[s2.index(sip1)])**2)/sigma
alpha1 = np.arctan((y1[s1.index(sip1)]-y2[s2.index(sip1)])/x1[s1.index(sip1)]-x2[s2.index(sip1)])*180./np.pi

d5 = np.sqrt((x1[s1.index(sip5b1)]-x2[s2.index(sip5b2)])**2+(y1[s1.index(sip5b1)]-y2[s2.index(sip5b2)])**2)/sigma
alpha5 = np.arctan((y1[s1.index(sip5b1)]-y2[s2.index(sip5b2)])/x1[s1.index(sip5b1)]-x2[s2.index(sip5b2)])*180./np.pi
print(angle[0],angle[2],d1,alpha1,d5,alpha5)
