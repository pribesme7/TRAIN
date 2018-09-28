import pickle as pk
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

with open('dicip1.pickle', 'rb') as handle:
    dicip1 = pk.load(handle)

with open('dicip5.pickle', 'rb') as handle:
    dicip5 = pk.load(handle)

print(dicip1["f"]["head"].diagonal(),len(dicip1["f"]["head"].diagonal()))
print(len(np.arange(-6.0, -6.0 + 0.1*len(dicip1["f"]["head"].diagonal()+0.1),0.1)),len(dicip1["f"]["head"].diagonal()))
#plt.plot(np.arange(-6.0, -6.0 + 0.1*len(dicip1["f"]["tail"].diagonal()+0.1),0.1),dicip1["f"]["head"].diagonal()-dicip1["f"]["tail"].diagonal())
#plt.show()
#print(dicip5)

step = 0.1
start = -6.0 - step
stop = 6.0 + step
    

onsep = np.arange(start + step,-step,step)
onsep2 = np.arange(-step,stop,step)
onsep = list(onsep) + list(onsep2)
onsep = [ '%.2f' % elem for elem in onsep ]
onsep = [ float(elem) for elem in onsep ]




ip1fhcosep5 = dicip1["f"]["head"][0,:]
ip1ftcosep5 = dicip1["f"]["tail"][0,:]
ip1fmcosep5 = dicip1["f"]["medium"][0,:]

ip1fhcosep1 = dicip1["f"]["head"][:,0]
ip1ftcosep1 = dicip1["f"]["tail"][:,0]
ip1fmcosep1 = dicip1["f"]["medium"][:,0]

ip1fhconsep15 = dicip1["f"]["head"]
ip1ftconsep15 = dicip1["f"]["tail"]
ip1fmconsep15 = dicip1["f"]["medium"]

plt.plot(onsep,ip1fhcosep1,label = 'head')
plt.plot(onsep,ip1ftcosep1,label = 'tail')
plt.plot(onsep,ip1fmcosep1,label = 'medium')
plt.legend()
plt.title('B1 IP1')
plt.xlabel('on_sep1')
plt.ylabel('closed orbit offset um')
plt.show()

plt.plot(onsep,ip1fhcosep5,label = 'head')
plt.plot(onsep,ip1ftcosep5,label = 'tail')
plt.plot(onsep,ip1fmcosep5,label = 'medium')
plt.legend()
plt.title('B1 IP1')
plt.xlabel('on_sep5')
plt.ylabel('closed orbit offset um')
plt.show()


ip1bhcosep5 = dicip1["b"]["head"][0,:]
ip1btcosep5 = dicip1["b"]["tail"][0,:]
ip1bmcosep5 = dicip1["b"]["medium"][0,:]

ip1bhcosep1 = dicip1["b"]["head"][:,0]
ip1btcosep1 = dicip1["b"]["tail"][:,0]
ip1bmcosep1 = dicip1["b"]["medium"][:,0]

ip1bhconsep15 = dicip1["b"]["head"]
ip1btconsep15 = dicip1["b"]["tail"]
ip1bmconsep15 = dicip1["b"]["medium"]

plt.plot(onsep,ip1bhcosep1,label = 'head')
plt.plot(onsep,ip1btcosep1,label = 'tail')
plt.plot(onsep,ip1bmcosep1,label = 'medium')
plt.legend()
plt.title('B2 IP1')
plt.xlabel('on_sep1')
plt.ylabel('closed orbit offset um')
plt.show()

plt.plot(onsep,ip1bhcosep5,label = 'head')
plt.plot(onsep,ip1btcosep5,label = 'tail')
plt.plot(onsep,ip1bmcosep5,label = 'medium')
plt.legend()
plt.title('B2 IP1')
plt.xlabel('on_sep5')
plt.ylabel('closed orbit offset um')
plt.show()

##########################################

ip5fhcosep5 = dicip5["f"]["head"][0,:]
ip5ftcosep5 = dicip5["f"]["tail"][0,:]
ip5fmcosep5 = dicip5["f"]["medium"][0,:]

ip5fhcosep1 = dicip5["f"]["head"][:,0]
ip5ftcosep1 = dicip5["f"]["tail"][:,0]
ip5fmcosep1 = dicip5["f"]["medium"][:,0]

ip5fhconsep15 = dicip5["f"]["head"]
ip5ftconsep15 = dicip5["f"]["tail"]
ip5fmconsep15 = dicip5["f"]["medium"]

plt.plot(onsep,ip5fhcosep1,label = 'head')
plt.plot(onsep,ip5ftcosep1,label = 'tail')
plt.plot(onsep,ip5fmcosep1,label = 'medium')
plt.legend()
plt.title('B1 IP5')
plt.xlabel('on_sep1')
plt.ylabel('closed orbit offset um')
plt.show()

plt.plot(onsep,ip5fhcosep5,label = 'head')
plt.plot(onsep,ip5ftcosep5,label = 'tail')
plt.plot(onsep,ip5fmcosep5,label = 'medium')
plt.legend()
plt.title('B1 IP5')
plt.xlabel('on_sep5')
plt.ylabel('closed orbit offset um')
plt.show()


ip5bhcosep5 = dicip5["b"]["head"][0,:]
ip5btcosep5 = dicip5["b"]["tail"][0,:]
ip5bmcosep5 = dicip5["b"]["medium"][0,:]

ip5bhcosep1 = dicip5["b"]["head"][:,0]
ip5btcosep1 = dicip5["b"]["tail"][:,0]
ip5bmcosep1 = dicip5["b"]["medium"][:,0]

ip5bhconsep15 = dicip5["b"]["head"]
ip5btconsep15 = dicip5["b"]["tail"]
ip5bmconsep15 = dicip5["b"]["medium"]

plt.plot(onsep,ip5bhcosep1,label = 'head')
plt.plot(onsep,ip5btcosep1,label = 'tail')
plt.plot(onsep,ip5bmcosep1,label = 'medium')
plt.legend()
plt.title('B2 IP5')
plt.xlabel('on_sep1')
plt.ylabel('closed orbit offset um')
plt.show()

plt.plot(onsep,ip5bhcosep5,label = 'head')
plt.plot(onsep,ip5btcosep5,label = 'tail')
plt.plot(onsep,ip5bmcosep5,label = 'medium')
plt.legend()
plt.title('B2 IP5')
plt.xlabel('on_sep5')
plt.ylabel('closed orbit offset um')
plt.show()



#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax = fig.gca(projection='3d')
#xx,yy = np.meshgrid(onsep,onsep)
#ax.scatter(xx,yy,ip1fhconsep15,label = 'head')
#ax.scatter(xx,yy,ip1ftconsep15,label = 'tail')
#ax.scatter(xx,yy,ip1fmconsep15,label = 'medium')

#ax.set_xlabel('on_sep1')
#ax.set_ylabel('on_sep5')
#ax.set_zlabel('closed orbit offset um')

#plt.legend()
#plt.show()"""



