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

filelist = [250,500,1000,1500,2000,3200,5000,10000,20000]
subscripts = ["","a"]
histogram = False
autocorr = False
poincare = True
attractor3d = False
points = False #legacy mode - look at photon counts
divs = True

if divs:
	def getkldiv(distr):
		#Kullback-Leibler
		distr = [float(x) / np.sum(distr) for x in distr] #Normalize
		uniformprob = np.ones(len(distr)) / len(distr)
		#kldiv = sum_i P_i log (P_i / Q_i)
		kldiv = 0
		for i in range(len(distr)):
			if distr[i] != 0:
				kldiv += distr[i] * np.log(distr[i] / uniformprob[i])
		print kldiv
		return kldiv

	def getksdiv(distr):
		distr = [float(x) / np.sum(distr) for x in distr] #Normalize
		distrcdf = np.cumsum(distr)
		uniformprob = np.ones(len(distr)) / len(distr)
		uniformcdf = np.cumsum(uniformprob)
		#ksdiv = sup_x | F_n(x) - F(x) |
		ksdiv = 0
		for i in range(len(distr)):
			ksdiv = max(ksdiv, abs(distrcdf[i] - uniformcdf[i]))
		print ksdiv
		return ksdiv

	kldivs = []
	ksdivs = []

for filename in filelist:
	print filename
	finvt = open(str(filename) + "vt.out","r")

	if histogram:
		finv = open(str(filename) + "v.out","r")
		voltages = []
		for line in finv:
			voltages.append(float(line))
		T = voltages.pop(0)
		print 'Voltages read!'

	pvoltages = [[],[]]
	T = float(finvt.readline())
	for line in finvt:
		vwp, vwpt = line.split()
		pvoltages[0].append(float(vwp))
		pvoltages[1].append(float(vwpt))

	window = int(Td / 4 * sampspersec)

	print 'Let the graphing begin!'
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
		tinterval = int(2 * transtime / samptime)
		vslice = voltages[:tinterval]
		iterlen = tinterval / 20
		mean = np.average(vslice)
		C = []
		for shift in range(iterlen):
			i = shift
			ans = 0.
			while i < tinterval:
				ans += (vslice[i - shift] - mean) * (vslice[i] - mean)
				i += 1
			ans /= (tinterval - shift)
			C.append(ans)

			#Progress
			if shift % (iterlen / 10) == 0:
				print int(shift * 100. / iterlen)
		C = [x / C[0] for x in C]
		tarr = [x * samptime for x in range(iterlen)]
		C.pop(0)
		tarr.pop(0)
		Cgraph = C[::-1] + [1] + C
		tgraph = [-x for x in tarr[::-1]] + [0] + tarr

		plt.figure(2)
		plt.title(str(filename))
		plt.plot(tgraph, Cgraph)
		print 'Autocorrelation done!'

	if poincare:
		pbinw = 0.008
		minp = 1.
		maxp = 5.
		ran = maxp - minp
		num = int(ran / pbinw)
		psec = [[0 for _ in range(num)] for x in range(num)]
		delay = Td / 4

		vertical = False #True if slope gets too high
		thickness = 3 #How many pbinws
		slicer = [[minp, maxp], [minp + 5. / 8 * ran, minp + 5. / 8 * ran]]
		if vertical:
			slopey = (slicer[0][1] - slicer[0][0]) / (slicer[1][1] - slicer[1][0])
			pslice = [0 for i in range(int((slicer[1][1] - slicer[1][0]) / pbinw))]
		else:
			slope = (slicer[1][1] - slicer[1][0]) / (slicer[0][1] - slicer[0][0])
			pslice = [0 for i in range(int((slicer[0][1] - slicer[0][0]) / pbinw))]

		for i in range(len(pvoltages[0])):
			vwp = pvoltages[0][i]
			vwpt = pvoltages[1][i]
			if vwp >= maxp or vwpt >= maxp or vwp <= minp or vwpt <= minp:
				continue
			psec[int((vwp - minp) / pbinw)][int((vwpt - minp) / pbinw)] += 1

			if vertical:
				orig = vwpt
				x = slicer[0][0] + slopey * (vwpt - slicer[1][0])
				if abs(int((x - minp) / pbinw) - int((vwp - minp) / pbinw)) <= thickness:
					pslice[int((vwpt - slicer[1][0]) / pbinw)] += 1
			else:
				orig = vwp
				y = slicer[1][0] + slope * (vwp - slicer[0][0])
				if abs(int((y - minp) / pbinw) - int((vwpt - minp) / pbinw)) <= thickness:
					pslice[int((vwp - slicer[0][0]) / pbinw)] += 1

		plt.figure(3)
		plt.subplot(121)
		plt.title(str(filename) + ", bin width = " + str(pbinw) + ", points = " + str(len(pvoltages[0])))
		plt.ylim([0,num])
		plt.xlim([0,num])
		plt.pcolormesh(np.transpose(np.array(psec)))
		scaledslicer = (np.array(slicer) - minp) / pbinw
		plt.plot(scaledslicer[0], scaledslicer[1], color = 'r')
		plt.subplot(122)
		plt.title("T = " + str(T) + ", thick = " + str(thickness))
		if not vertical:
			plt.xlim([0,len(pslice)])
			plt.bar(range(len(pslice)), pslice)
		else:
			plt.ylim([0,len(pslice)])
			plt.barh(range(len(pslice)), pslice)			
		print 'Poincare section done!'

		if divs:
			trimmed = np.trim_zeros(pslice)
			kldivs.append(getkldiv(trimmed))
			ksdivs.append(getksdiv(trimmed))

	if attractor3d:
		bigplot = [[],[],[]]
		for i in range(2 * window, len(voltages) / 100):
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

	#plt.show()

if divs:
	plt.figure(6)
	plt.subplot(211)
	plt.title("Kullback-Leibler distance")
	plt.plot(kldivs)
	plt.subplot(212)
	plt.title("Kolmogorov-Smirnov distance")
	plt.plot(ksdivs)
	plt.show()

print 'Program done'
