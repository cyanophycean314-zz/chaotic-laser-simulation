#Graph data of interarrival times

import sys
import matplotlib.pyplot as plt
import numpy as np

#The (basically unchangeable conditions)
dt = 0.00002 #Time interval (around the size of mu for lambda0 * Td = 3200)
samptime = 0.00002 #How often to take a sample
sampspersec = 1 / samptime #Inverse
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
transtime = 0.1
phi = np.pi / 4 #Filter phase displacement
veps = 0.07
psliceval = 3.1415926535
lengthofgraph = 10000

#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.

filename = ''.join(sys.argv[1:])

fin = open(filename, 'r')
times = []
for line in fin:
	times.append(float(line))
actualtime = np.cumsum(times)/151.515e6

T = actualtime[-1]
binw = Td / 4
w = Td / 4

def getmovingwindow(photonpops):
	binno = int(T / binw) + 1
	print binno
	movingwindow = [0] * (binno)
	#print binno

	lowcnt = 0 #binsearch(photonpops, transcount * Td * 0.99)
	highcnt = 0 #binsearch(photonpops, transcount * Td)
	for i in range(binno):
		#Sliding count - linear, suggested by Joe
		while lowcnt < len(photonpops) and photonpops[lowcnt] < i * binw:
			lowcnt += 1
		while highcnt < len(photonpops) and photonpops[highcnt] < i * binw + w:
			highcnt += 1
		movingwindow[i] = highcnt - lowcnt

	#print movingwindow
	#Toss out initial transient phase
	print 'Histogram complete!'
	timegraph = [x * binw for x in range(binno)]

	return timegraph, movingwindow

tg, mw = getmovingwindow(actualtime)

plt.figure(1)
plt.subplot(211)
plt.title('lambda0 = 2.76MHz, lambda0Td = 4785')
plt.xlim([0,T])
plt.xlabel("Time (s)")
plt.ylabel("Photons in moving window")
plt.plot(tg, mw)
plt.subplot(212)
plt.title("Histogram of above graph")
plt.xlabel("Photons in moving window")
plt.ylabel("Bin count")
plt.hist(mw, bins = 30)


pbinw = 4
maxp = 800
num = int(maxp / pbinw)
psec = [[0 for x in range(num)] for y in range(num)]
for ptime in poincaretimes:
	nwp = movingwindowcalc(ptime)
	nwpt = movingwindowcalc(ptime - Td / 4)
	if nwp >= maxp or nwpt >= maxp:
		continue
	#print str(nwp) + "," + str(nwpt) + "," + str(maxp)
	psec[int(nwp / pbinw)][int(nwpt / pbinw)] += 1

print 'Point-care section done! Ha get it?'
plt.figure(5)
plt.pcolor(np.array(psec))

plt.show()
