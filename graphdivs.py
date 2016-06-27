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

filelist = [10000]#[250,500,1000,1500,2000,3200,5000,10000]#,10000,20000]
subscripts = list("abcdefghij")
#filelist = ["detBIG"]
#subscripts = [""]
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
	pbinw = 0.006
	minp = 1.
	maxp = 5.
	ran = maxp - minp
	num = int(ran / pbinw)
	psec = [[0 for _ in range(num)] for x in range(num)]
	delay = Td / 4

	vertical = False #True if slope gets too high
	thickness = 2 #How many pbinws
	slicer = [[minp, maxp], [maxp, minp]]
	if vertical:
		slopey = (slicer[0][1] - slicer[0][0]) / (slicer[1][1] - slicer[1][0])
		pslice = [0 for i in range(int((slicer[1][1] - slicer[1][0]) / pbinw))]
	else:
		slope = (slicer[1][1] - slicer[1][0]) / (slicer[0][1] - slicer[0][0])
		pslice = [0 for i in range(int((slicer[0][1] - slicer[0][0]) / pbinw))]
	for lett in subscripts:
		print str(filename) + lett
		finvt = open(str(filename) + lett + "vt.out","r")

		pvoltages = [[],[]]
		T = float(finvt.readline())
		for line in finvt:
			vwp, vwpt = line.split()
			pvoltages[0].append(float(vwp))
			pvoltages[1].append(float(vwpt))

		window = int(Td / 4 * sampspersec)
		print 'Let the graphing begin!'

		for i in range(len(pvoltages[0])):
			vwp = pvoltages[0][i]
			vwpt = pvoltages[1][i]
			if int((vwp - minp) / pbinw) >= num or int((vwpt - minp) / pbinw) >= num or vwp <= minp or vwpt <= minp:
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

	plt.figure(3, figsize = (25,10))
	plt.subplot(121)
	plt.title(str(filename) + ", bin width = " + str(pbinw) + ", points = " + str(np.sum(psec)))
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
	plt.savefig(str(filename) + ".png")

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
	#plt.show()

if divs:
	plt.figure(6)
	plt.subplot(211)
	plt.title("Kullback-Leibler distance")
	plt.plot(kldivs)
	plt.subplot(212)
	plt.title("Kolmogorov-Smirnov distance")
	plt.plot(ksdivs)
	plt.savefig("divs.png")

print 'Program done'
