#Graph the files

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time

#The (basically unchangeable conditions)
dt = 0.000005 #Time interval (around the size of mu for lambda0 * Td = 3200)
samptime = 0.00002 #How often to take a sample
sampspersec = 1 / samptime #Inverse
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
transtime = 0.1
phi = np.pi / 4 #Filter phase displacement

#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
T = 5 #seconds to simulate
deterministic = False

filelist = [3200]
histogram = False
autocorr = False
poincare = True
attractor3d = False
points = False #legacy mode - look at photon counts

for filename in filelist:
	finv = open(str(filename) + "v.out","r")
	finx = open(str(filename) + "xs.out","r")

	voltages = []
	for line in finv:
		voltages.append(float(line))
	T = voltages.pop(0)

	window = int(Td / 4 * sampspersec)
	smoothedvoltages = [np.mean(voltages[i - window: i]) for i in range(window, len(voltages))]
	rawvoltages = voltages
	voltages = smoothedvoltages
	timegraph = [transtime + Td / 4 + x * samptime for x in range(len(smoothedvoltages))]

	poincaretimes = []
	for line in finx:
		poincaretimes.append(float(line))
		#print rawvoltages[int((float(line) - transtime) * sampspersec)]
		#print voltages[int((float(line) - transtime - Td / 4) * sampspersec)]
	while max(poincaretimes) > timegraph[len(timegraph) - 1]:
		poincaretimes.pop()

	#Graph the stuff
	if histogram:
		plt.figure(1)
		plt.subplot(211)
		plt.title(str(filename))
		plt.xlim([0, T - transtime])
		plt.plot(timegraph, voltages)
		plt.subplot(212)
		plt.hist(voltages, bins = 100)
		print 'Histogram done!'

	if autocorr:
		tinterval = 2 * transtime
		vslice = voltages[:tinterval]
		iterlen = tinterval / 20
		mean = np.average(vslice)
		C = []
		for shift in range(iterlen):
			i = shift
			ans = 0.
			while i < tinterval:
				ans += (vslice[i - shift] - mean) * (vslice[j] - mean)
				i += 1
			ans /= (tinterval - shift)
			C.append(ans)

			#Progress
			if shift % (iterlen / 10) == 0:
				print int(shift * 100. / iterlen)
		C = [x / C[0] for x in C]
		tarr = np.array(range(iterlen)) * samptime
		C.pop(0)
		tarr.pop(0)
		Cgraph = C[::-1] + [1] + C
		tgraph = tarr[::-1] + [0] + tarr

		plt.figure(2)
		plt.title(str(filename))
		plt.plot(tgraph, Cgraph)
		print 'Autocorrelation done!'

	if poincare:
		pbinw = 0.05
		minp = 0.5
		maxp = 6.
		ran = maxp - minp
		num = int(ran / pbinw)
		psec = [[0 for _ in range(num)] for x in range(num)]
		delay = Td / 4

		slicer = [[minp + ran / 4, minp + 3 * ran / 4], [minp, maxp]]
		slope = (slicer[1][1] - slicer[1][0]) / (slicer[0][1] - slicer[0][0])
		pslice = [0 for i in range(int((slicer[0][1] - slicer[0][0]) / pbinw))]
		for ptime in poincaretimes:
			if ptime > timegraph[0] + delay:
				pcorr = ptime - timegraph[0] - Td / 4
				vwp = voltages[int(pcorr * sampspersec)]
				vwpt = voltages[int((pcorr - delay) * sampspersec)]

				if vwp >= maxp:
					vwp -= pbinw / 2
				elif vwp <= minp:
					vwp += pbinw / 2
				if vwpt >= maxp:
					vwpt -= pbinw / 2
				elif vwpt <= minp:
					vwpt += pbinw / 2

				psec[int((vwp - minp) / pbinw)][int((vwpt - minp) / pbinw)] += 1

				#1 bin margin of error
				y = slicer[1][0] + slope * (vwp - slicer[0][0])
				#print "{} {}".format(int(y / pbinw), int(vwpt / pbinw))
				if abs(int(y / pbinw) - int(vwpt / pbinw)) <= 1:
					pslice[int((vwp - slicer[0][0]) / pbinw)] += 1

		plt.figure(3)
		plt.subplot(211)
		plt.title(str(filename) + ", bin width = " + str(pbinw))
		plt.ylim([0,num])
		plt.xlim([0,num])

		plt.pcolor(np.array(psec))
		scaledslicer = (np.array(slicer) - minp) / pbinw
		plt.plot(scaledslicer[0], scaledslicer[1], color = 'r')
		plt.subplot(212)
		plt.scatter(range(len(pslice)), pslice)

		print 'Poincare section and slice done!'

	if attractor3d:
		bigplot = [[],[],[]]
		for i in range(2 * window, len(voltages)):
			bigplot[0].append(voltages[i])
			bigplot[1].append(voltages[i - window])
			bigplot[2].append(voltages[i - 2 * window])

		fig = plt.figure(4)
		ax = fig.add_subplot(111, projection = '3d')
		ax.plot(bigplot[0], bigplot[1], bigplot[2])

		print '3D attractor done!'

	if points:
		finpop = open(str(filename) + "pop.out","r")
		photonpops = []
		for line in finpop:
			photonpops.append(float(line))

		def binsearch(lst, key):
			low = 0
			top = len(lst) - 1
			
			while low < top:
				#Loop exits when low = top
				med = (low + top) / 2
				if key < lst[med]:
					top = med
				else:
					low = med + 1
			return low

		def movingwindowcalc(t):
			#print binno
			lowerbound = binsearch(photonpops, t - Td / 4)
			upperbound = binsearch(photonpops, t)
			return upperbound - lowerbound

		pbinw = 4
		maxp = 800
		num = maxp / pbinw
		psec = [[0 for x in range(num)] for y in range(num)]
		for ptime in poincaretimes:
			if ptime >= transtime + Td / 4:
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

	plt.show()

print 'Program done'
