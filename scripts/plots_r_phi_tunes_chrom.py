import numpy as np
import matplotlib.pyplot as plt
from pylab import rcParams
import pylab
rcParams['figure.figsize'] = 20, 5

#path = 'RESULTS/MDcycle/collision_B_IP15/fort.ps.hoff_bIP1'
path2 = 'RESULTS/MDcycle_new_no_oct_no_disp/collision_41_std_IP1258/tune_f.list'
path7 = 'RESULTS/MDcycle_new_no_oct_no_disp/collision_41_std_IP1258/tune_b.list'
path13 = 'RESULTS/MDcycle/provA/tune_f.list'


#path3 = 'RESULTS/MDcycle/collision_B_IP15/hoff_f.IP1'
#path4 = 'RESULTS/MDcycle/collision_B_IP15/hoff_b.IP1'
#path5 = 'RESULTS/BCMS_alternative.in/hoff_f.IP1'
#path6 = 'RESULTS/BCMS_alternative.in/hoff_b.IP1'
#slotnb, x,xp,r,phi = np.loadtxt(path,unpack = True,skiprows=0, usecols=[0, 1,2,3,4])
slot,bc,th,tv,ch,cv = np.loadtxt(path2,unpack = True,skiprows=0, usecols=[0, 1,2,3,4,5])
#slot3,hb1 = np.loadtxt(path3,unpack = True,skiprows=0, usecols=[0, 1])
slot,bc2,th2,tv2,ch2,cv2 = np.loadtxt(path7,unpack = True,skiprows=0, usecols=[0, 1,2,3,4,5])
#slot4,hb2 = np.loadtxt(path4,unpack = True,skiprows=0, usecols=[0, 1])
#slot5,hb5 = np.loadtxt(path5,unpack = True,skiprows=0, usecols=[0, 1])
#slot6,hb6 = np.loadtxt(path6,unpack = True,skiprows=0, usecols=[0, 1])
listpacman8 = [808,809,810,811,812,813,814,815,816,817
,818,819,820,821,822,823,824,825,826,827
,828,829,830,831,832,833,834,835,836,837
,838,839,840,841,842,843,844,845,846,847
,848,849,850,851,852,853,854,855,856,857
,858,859,860,861,862,863,864,865,866,867
,868,869,870,871,872,873,874,875,876,877
,878,879,887,888,889,890,891,892,893,906
,907,908,909,910,911,912,913,914,915,916
,917,918,919,920,921,922,923,924,925,926
,927,928,929,930,931,932,933,934,935,936
,937,938,939,940,941,942,943,944,945,946
,947,948,990,991,992,993,994,995,996,997
,998,999,1000,1001,1074,1075,1076,1077,1078,1079
,1080,1153,1154,1155,1156,1157,1158,1159,1232,1233
,1234,1235,1236,1237,1238,1330,1331,1332,1333,1334
,1335,1336,1337,1338,1339,1340,1341,1414,1415,1416
,1417,1418,1419,1420,1493,1494,1495,1496,1497,1498
,1499,1572,1573,1574,1575,1576,1577,1578]

listpacman2 = [171,172,173,174,175,176,177,250,251,252
,253,254,255,256,329,330,331,332,333,334
,335,408,409,410,411,412,413,414,415,416
,511,512,513,514,515,516,517,590,591,592
,593,594,595,596,669,670,671,672,673,674
,675,748,749,750,751,752,753,754,755,756
,808,809,810,887,888,889,990,991,992,1069
,1070,1071,1148,1149,1150,1227,1228,1229,1330,1331
,1332,1409,1410,1411,1488,1489,1490,1567,1568,1569
,1702,1703,1704,1781,1782,1783,1884,1885,1886,1963
,1964,1965,2042,2043,2044,2121,2122,2123,2224,2225
,2226,2303,2304,2305,2382,2383,2384,2461,2462,2463
,2596,2597,2598,2599,2600,2601,2602,2603,2604,2605
,2606,2607,2608,2609,2610,2611,2612,2613,2614,2615
,2616,2617,2618,2619,2620,2621,2622,2623,2624,2625
,2626,2627,2628,2629,2630,2631,2632,2633,2634,2635
,2636,2637,2638,2639,2640,2641,2642,2643,2644,2645
,2646,2647,2648,2649,2650,2651,2652,2653,2654,2655
,2656,2657,2658,2659,2660,2661,2662,2663,2664,2665
,2666,2667,2685,2686,2687,2688,2689,2690,2691,2692
,2693,2694,2695,2696,2697,2698,2699,2700,2701,2702
,2703,2704,2705,2706,2707,2708,2709,2710,2711,2712
,2713,2714,2715,2716,2717,2718,2719,2720,2721,2722
,2723,2724,2725,2726,2727,2740,2741,2742,2743,2744
,2745,2746,2778,2779,2780,2857,2858,2859,2936,2937
,2938,3015,3016,3017,3118,3119,3120,3197,3198,3199
,3276,3277,3278,3355,3356,3357 ]



"""
plt.plot(slot3,hb1,'ro',label = 'B1',alpha = 0.5)
plt.plot(slot4,hb2,'bo',label='B2',alpha = 0.5)
plt.legend()
plt.xlabel('slot id')
plt.ylabel('Horizontal orbit distortions (um)')
plt.title('IP1')
plt.show()

plt.plot(slot5,hb6,'ro',label = 'B1',alpha = 0.5)
plt.plot(slot5,hb5,'bo',label='B2',alpha = 0.5)
plt.legend()
plt.xlabel('slot id')
plt.ylabel('Horizontal orbit distortions (um)')
plt.title('IP1 (BCMS Alternative)')
plt.show()
"""


plt.plot(slot[25:],th[25:],'bo',label='B1', alpha = 0.5)
plt.plot(slot[25:],th2[25:],'ro',label='B2',alpha = 0.5)
#plt.title('Forward Beam, IP1, std, Horizontal')
plt.xlabel('bunch id',fontsize=20)
plt.ylabel('Horizontal Tune',fontsize=20)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.legend()
plt.show()


plt.plot(slot[25:],tv[25:],'bo',label='B1', alpha = 0.5)
plt.plot(slot[25:],tv2[25:],'ro',label='B2',alpha = 0.5)
#plt.title('Forward Beam, IP1, std, Vertical')
plt.xlabel('bunch id',fontsize=20)
plt.ylabel('Vertical Tune',fontsize=20)
plt.legend()
plt.grid(True)
plt.xticks(fontsize=18, rotation=0)
plt.yticks(fontsize=18, rotation=0)
plt.show()

plt.plot(slot[25:],ch[25:],'ro',label='B1',alpha = 0.5)
plt.plot(slot[25:],ch2[25:],'bo',label='B2',alpha = 0.5)
#plt.title('Forward Beam, IP1, std')
plt.xlabel('bunch id',fontsize=15)
plt.ylabel('Horizontal Chromaticity',fontsize=15)
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.legend()
plt.grid(True)
plt.show()

plt.plot(slot[25:],cv[25:],'ro',label='B1',alpha = 0.5)
plt.plot(slot[25:],cv2[25:],'bo',label='B2',alpha = 0.5)
#plt.title('Forward Beam, IP1, std')
plt.xlabel('bunch id',fontsize=15)
plt.ylabel('Vertical Chromaticity',fontsize=15)
plt.xticks(fontsize=15, rotation=0)
plt.yticks(fontsize=15, rotation=0)
plt.legend()
plt.show()
"""
fig, ax1 = plt.subplots()

ax1.plot(slotnb, r, 'b-')
ax1.set_xlabel('bunch id')
# Make the y-axis label, ticks and tick labels match the line color.
ax1.set_ylabel('r (sigma)', color='b')
ax1.tick_params('y', colors='b')

ax2 = ax1.twinx()
ax2.plot(slotnb, phi, 'r.')
ax2.set_ylabel('phi', color='r')
ax2.tick_params('y', colors='r')

fig.tight_layout()
plt.show()"""

"""

plt.plot(slot[25:],th[25:],'bo',label='B1', alpha = 0.5)
#plt.plot(slot[25:],th2[25:],'ro',label='B2',alpha = 0.5)
for i in listpacman8:

	plt.plot(i,th[list(slot).index(i)],color = 'green',marker='o',alpha = 0.5)

for i in listpacman2:
	plt.plot(i,th[list(slot).index(i)],color = 'pink',marker = 'o',alpha = 0.5)
plt.xlabel('bunch id')
plt.ylabel('Horizontal tunes')
plt.legend()
plt.show()

plt.plot(slot[25:],tv[25:],'bo',label='B1', alpha = 0.5)
#plt.plot(slot[25:],tv2[25:],'ro',label='B2',alpha = 0.5)
for i in listpacman8:
	
	plt.plot(i,tv[list(slot).index(i)],color = 'green',marker='o',alpha = 0.5)

for i in listpacman2:
	plt.plot(i,tv[list(slot).index(i)],color = 'pink',marker = 'o',alpha = 0.5)

plt.xlabel('bunch id')
plt.ylabel('Vertical tunes')
plt.legend()
plt.show()

plt.plot(slot[25:],ch[25:],'bo',label='B1',alpha = 0.5)
#plt.plot(slot[25:],ch2[25:],'ro',label='B2',alpha = 0.5)

for i in listpacman8:
	
	plt.plot(i,ch[list(slot).index(i)],color = 'green',marker='o',alpha = 0.5)

for i in listpacman2:
	plt.plot(i,ch[list(slot).index(i)],color = 'pink',marker = 'o',alpha = 0.5)
plt.xlabel('bunch id')
plt.ylabel('Horizontal chromaticities')
plt.legend()
plt.show()

plt.plot(slot[25:],cv[25:],'bo',label='B1',alpha = 0.5)
#plt.plot(slot[25:],cv2[25:],'ro',label='B2',alpha = 0.5)
for i in listpacman8:
	plt.plot(i,cv[list(slot).index(i)],color = 'green',marker='o',alpha = 0.5)

for i in listpacman2:
	plt.plot(i,cv[list(slot).index(i)],color = 'pink',marker = 'o',alpha = 0.5)
plt.xlabel('bunch id')
plt.ylabel('Vertical chromaticities')
plt.legend()
plt.show()"""
