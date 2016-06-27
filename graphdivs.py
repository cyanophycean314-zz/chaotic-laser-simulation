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

filelist = [10,100,250,500,1000,1500,2000,3200,5000,10000,20000,30000]
subs = list("abcdefghijklm")
#filelist = ["detx" + _ for _ in ["005","008","01","03","05","07","1"]]
#subs = [""]#["BIG"] + list("abcdefgh")
deterministic = True

histogram = False
autocorr = False
poincare = True
attractor3d = False
points = False #legacy mode - look at photon counts
divs = True

if divs:
	def getkldiv(distr, otherdist = "uniform"):
		#Kullback-Leibler
		if otherdist == "uniform":
			distr = np.trim_zeros(distr)
			compprob = np.ones(len(distr)) / len(distr)
		else:
			compprob = [float(x) / np.sum(otherdist) for x in otherdist]
		distr = [float(x) / np.sum(distr) for x in distr] #Normalize
		#kldiv = sum_i P_i log (P_i / Q_i)
		kldiv = 0
		for i in range(len(distr)):
			if distr[i] * compprob[i] != 0:
				kldiv += distr[i] * np.log(distr[i] / compprob[i])
		print kldiv
		return kldiv

	def getksdiv(distr, otherdist = "uniform"):
		if otherdist == "uniform":
			distr = np.trim_zeros(distr)
			compprob = np.ones(len(distr)) / len(distr)
			compcdf = np.cumsum(compprob)
		else:
			compprob = [float(x) / np.sum(otherdist) for x in otherdist]
			compcdf = np.cumsum(compprob)
		distr = [float(x) / np.sum(distr) for x in distr] #Normalize
		distrcdf = np.cumsum(distr)
		#ksdiv = sup_x | F_n(x) - F(x) |
		ksdiv = 0
		for i in range(len(distr)):
			ksdiv = max(ksdiv, abs(distrcdf[i] - compcdf[i]))
		print ksdiv
		return ksdiv

	kldivs = []
	ksdivs = []

#Calculate pure deterministic version
pbinw = 0.01
minp = 1.
maxp = 5.
ran = maxp - minp
num = int(ran / pbinw)
delay = Td / 4
vertical = False #True if slope gets too high
thickness = 2 #How many pbinws
slicer = [[minp + 2 * ran / 5, minp + 3 * ran / 5], [maxp, minp]]
if vertical:
	slopey = (slicer[0][1] - slicer[0][0]) / (slicer[1][1] - slicer[1][0])
else:
	slope = (slicer[1][1] - slicer[1][0]) / (slicer[0][1] - slicer[0][0])

noises = []
for fileno in range(len(filelist) + 1):
	if fileno == 0:
		filename = "det"
		subscripts = ["BIG"] + list("abcdefgh")
	else:
		filename = filelist[fileno - 1]
		subscripts = subs

	#Assign noise
	if deterministic:
		if filename[:4] == "detx":
			noise = float("0." + filename[4:])
		else:
			noise = 0
	else:
		noise = 1 / np.sqrt(float(filename))

	if fileno != 0:
		noises.append(noise)

	#Reset data
	psec = [[0 for _ in range(num)] for x in range(num)]
	if vertical:
		pslice = [0 for i in range(int((slicer[1][1] - slicer[1][0]) / pbinw))]
	else:
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
		#print 'Let the graphing begin!'

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

	plt.figure(10+fileno, figsize = (25,10))
	plt.subplot(121)
	plt.title(str(filename) + ", bin width = " + str(pbinw) + ", points = " + str(np.sum(psec)))
	plt.ylim([0,num])
	plt.xlim([0,num])
	plt.pcolormesh(np.transpose(np.array(psec)))
	scaledslicer = (np.array(slicer) - minp) / pbinw
	plt.plot(scaledslicer[0], scaledslicer[1], color = 'r')
	plt.subplot(122)
	plt.title("T = " + str(T) + ", thick = " + str(thickness) + ", noise = " + str(noise))
	if not vertical:
		plt.xlim([0,len(pslice)])
		plt.bar(range(len(pslice)), pslice)
	else:
		plt.ylim([0,len(pslice)])
		plt.barh(range(len(pslice)), pslice)			
	print 'Poincare section done!'
	plt.savefig(str(filename) + ".png")

	if fileno == 0:
		pslicedet = pslice

	if divs and fileno != 0:
		kldivs.append(getkldiv(pslice, pslicedet))
		ksdivs.append(getksdiv(pslice, pslicedet))

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
	plt.plot(noises, kldivs)
	plt.subplot(212)
	plt.title("Kolmogorov-Smirnov distance")
	plt.plot(noises, ksdivs)
	plt.savefig("divs.png")

	fout = open("divs.out","w")
	for i in range(len(noises)):
		fout.write("{:6f} {:6f} {:6f}".format(noises[i], kldivs[i], ksdivs[i]))

print 'Program done'
