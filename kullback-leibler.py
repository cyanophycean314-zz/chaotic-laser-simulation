#Kullback Leibler
#Finds and graphs the Kullbach Liebler divergence for different lambda0

'''The Plan
	Read in the data
	Construct the pdf (with appropriate number of bins)
	Apply the Kullback-Leibler Metric
'''

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time

#Conditions
dt = 0.0000005 #Time interval (around the size of mu for lambda0 * Td = 3200)
Td = 0.001734 #Time delay
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
lambda0timesTd = 3200 #Metric given in Aaron's paper
lambda0 = lambda0timesTd / Td
mu = 1 / lambda0 #Poisson interarrival time average
phi = np.pi / 4 #Filter phase displacement
transcount = 100 #How many Td's to wait for transients to die
n = 150000000 #Photons to generate
xeps = 0.003 #X epsilon
eps = 0.00001 #other epsilon

filelist = ["12,5", 25, 50, 100, 200, 400, 600, 800, 1200, 1600, 2400, 3200]
binnings =  [1,1,1,1,2,2,2,2,2,2,2,2]

#Compare with the uniform
######################################################
def getkldiv(mw, dFU, binwidth):
	#Kullback-Leibler
	mwhistprob = np.histogram(movingwindow, bins = range(0, max(mw) + binwidth, binwidth))[0] / float(len(movingwindow))
	#Make the uniform
	uniformprob = np.ones(len(mwhistprob)) / len(mwhistprob)
	#kldiv = sum_i P_i log (P_i / Q_i)
	binnumber = len(range(0, max(mw), binwidth))
	kldiv = 0
	for i in range(binnumber):
		if dFU:
			kldiv += uniformprob[i] * np.log(uniformprob[i] / mwhistprob[i])
		else:
			kldiv += mwhistprob[i] * np.log(mwhistprob[i] / uniformprob[i])
	print kldiv
	'''
	if dFU:
		plt.figure(1)
		plt.subplot(211)
		plt.title('lambda0timesTd = ' + str(lambda0timesTd))
		plt.bar(range(binnumber),mwhistprob)
		plt.subplot(212)
		plt.bar(range(binnumber),uniformprob)
		plt.show()
	'''
	
	return kldiv	

def getksdiv(mw):
	mw = movingwindow
	probno = max(mw)
	mwhistprob = np.histogram(mw, bins = range(probno + 1))[0] / float(len(mw))
	mwhistcdf = np.cumsum(mwhistprob)
	#Kolmogorov-Smirnov statistic
	uniformprob = np.ones(probno) / probno
	uniformcdf = np.cumsum(uniformprob)
	#ksdiv = sup_x | F_n(x) - F(x) |
	ksdiv = 0
	for i in range(probno):
		ksdiv = max(ksdiv, abs(mwhistcdf[i] - uniformcdf[i]))
	print ksdiv
	'''
	plt.figure(1)
	plt.subplot(211)
	plt.title('n = ' + str(n))
	plt.plot(mwhistcdf)
	plt.subplot(212)
	plt.plot(uniformcdf)
	plt.show()
	'''
	return ksdiv

klgraph1 = []
klgraph2 = []
ksgraph = []
for fileno in range(len(filelist)):
	simulation = filelist[fileno]

	#Read the files
	######################################################
	finpop = open(str(simulation) + "pop.out", 'r')
	finxs = open(str(simulation) + "xs.out", 'r')
	photonpops = []
	poincaretimes = []
	for line in finpop:
		photonpops.append(line)
	photonpops.pop() #Last line is a bit buggy usually
	runinfo = photonpops.pop(0).split() #First line is the info section
	t = 0
	for line in finxs:
		poincaretimes.append(float(line))
	for i in range(len(photonpops)):
		photonpops[i] = float(photonpops[i])
	n = int(float(runinfo[2]))
	lambda0timesTd = float(runinfo[1])
	T = float(runinfo[0])

	print 'Done reading file = ' + str(simulation) + '!'

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

	#Generate the histogram
	######################################################
	w = Td / 4
	binw = w / 40
	offset = int(transcount * Td / binw)
	binno = int(T / binw) - offset
	movingwindow = [0] * binno
	divFromUniform = True #P = uniform dist

	lowcnt = 0 #binsearch(photonpops, transcount * Td * 0.99)
	highcnt = 0 #binsearch(photonpops, transcount * Td)
	for i in range(offset, offset + binno):
		#Sliding count - linear, suggested by Joe
		while lowcnt < len(photonpops) and photonpops[lowcnt] < i * binw:
			lowcnt += 1
		while highcnt < len(photonpops) and photonpops[highcnt] < i * binw + w:
			highcnt += 1
		movingwindow[i - offset] = highcnt - lowcnt

	klgraph1.append(getkldiv(movingwindow, divFromUniform, binnings[fileno]))
	klgraph2.append(getkldiv(movingwindow, not divFromUniform, binnings[fileno]))
	ksgraph.append(getksdiv(movingwindow))

filelist[0] = 12.5
plt.figure(2)

plt.subplot(311)
plt.title("Kullback-Leibler Divergence of Pdf from Uniform")
plt.ylabel("Arbitrary Units (a.u.)")
plt.plot(filelist, klgraph1)

plt.subplot(312)
plt.title("Kullback-Leibler Divergence of Uniform from Pdf")
plt.ylabel("Arbitrary Units (a.u.)")
plt.plot(filelist, klgraph2)

plt.subplot(313)
plt.title("Kolmogorov-Smirnov Statistic")
plt.xlabel("Photon rate * Td")
plt.ylabel("Arbitrary Units (a.u.)")
plt.plot(filelist, ksgraph)

plt.show()