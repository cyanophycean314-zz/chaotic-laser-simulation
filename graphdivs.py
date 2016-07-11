#Graph the files

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time

#The (basically unchangeable conditions)
dt = 0.00002 #Time interval (around the size of mu for lambda0 * Td = 3200)
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

#[100,250,500,1000,2000,3200,5000,10000,20000,30000]
filelist = ["5000_" + str(_ * 10)  for _ in range(10)]
subs = ['']
subsubscripts = [''] #+ list("ab")
deterministic = False

autocorr = False
poincare = True
attractor3d = False
points = False #legacy mode - look at photon counts
divs = False

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
pbinw = 0.005
minp = 1.
maxp = 5.
ran = maxp - minp
num = int(ran / pbinw)
delay = Td / 4
vertical = False #True if slope gets too high
thickness = 1 #How many pbinws
#slicer = [[minp + 1.3 * ran / 8, minp + 7.3 * ran / 8], [minp, maxp]]
slicer = [[minp, maxp], [minp + 5.05 * ran / 8, minp + 6 * ran / 8]]
if vertical:
	slopey = (slicer[0][1] - slicer[0][0]) / (slicer[1][1] - slicer[1][0])
else:
	slope = (slicer[1][1] - slicer[1][0]) / (slicer[0][1] - slicer[0][0])

noises = []
for fileno in range(len(filelist) + 1):
	if fileno == 0:
		if not divs:
			continue
		filename = "det"
		subscripts = ['com']
	else:
		filename = filelist[fileno - 1]
		subscripts = subs

	#Assign noise
	if divs:
		if deterministic or fileno == 0:
			if len(filename) >= 4 and filename[:4] == "detx":
				noise = float("0." + filename[4:])
			elif len(filename) >= 2 and filename[-2:] == "NT":
				noise = noise = 1 / np.sqrt(float(filename[:-2]))
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
	if divs:
		plt.title("T = " + str(T) + ", thick = " + str(thickness) + ", noise = " + str(noise))
	else:
		plt.title("T = " + str(T) + ", thick = " + str(thickness))
	if not vertical:
		plt.xlim([len(pslice) / 2,len(pslice)])
		plt.bar(range(len(pslice)), pslice)
	else:
		plt.ylim([len(pslice) / 2,len(pslice)])
		plt.barh(range(len(pslice)), pslice)			
	print 'Poincare section done!'
	plt.savefig(str(filename) + ".png")

	if fileno == 0:
		pslicedet = pslice

	if divs and fileno != 0:
		kldivs.append(getkldiv(pslice[len(pslice)/2:], pslicedet[len(pslice)/2:]))
		ksdivs.append(getksdiv(pslice[len(pslice)/2:], pslicedet[len(pslice)/2:]))

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
		fout.write("{:6f} {:6f} {:6f}\n".format(noises[i], kldivs[i], ksdivs[i]))

print 'Program done'
