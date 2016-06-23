#!/usr/bin/env python
# Simulate the differential equations in Hagerstrom (2015)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import random
import time
import math


'''
Outline

Generate photon times
Loop
	Get a photon
	Go through modulator
	Go through filter
'''

#The (basically unchangeable conditions)
dt = 0.0000005 #Time interval (around the size of mu for lambda0 * Td = 3200)
Td = 0.00173431292 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
phi = np.pi / 4 #Filter phase displacement
xeps = 0.003 #X epsilon
eps = 0.00001 #other epsilon
#Simulation parameters
betatimesTd = 8.86773684211 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
lambda0timesTd = 3200 #Metric given in Aaron's paper
lambda0 = lambda0timesTd / Td
mu = 1 / lambda0 #Poisson interarrival time average
transcount = 100 #How many Td's to wait for transients to die
n = 1000000 #Photons to generate
#Histogram parameters
w = Td / 4
binw = w / 20
offset = int(transcount * Td / binw)
maxp = 1
pbinw = 1

deterministic = True
filelist = ["lambda0Td=" + str(x) for x in [1000,1500,2000,2500,3000,3500,4000,5000,10000]]
#filelist = ["lambda0Td=10000"]
histogram = False
autocorr = False
poincare = True
attractor3d = False
variance = False

if not deterministic:
	for lambda0timesTd in [8000]:
		print lambda0timesTd
		lambda0 = lambda0timesTd / Td
		mu = 1 / lambda0 #Poisson interarrival time average
		n = 10000 * lambda0timesTd

		foutpop = open('lambda0Td=' + str(lambda0timesTd) + 'pop.out','w')
		foutx = open('lambda0Td=' + str(lambda0timesTd) + 'xs.out','w') #x-simplified
		foutv = open('lambda0Td=' + str(lambda0timesTd) + 'v.out','w') #x-simplified
		#Generate photon times!
		taus = np.random.exponential (mu , n) #Interarrival times
		photontimes = np.cumsum(taus) #Cumulative sum, times of arrival at modulator
		photonpops = [] #Times of photon arivals at counter!
		T = photontimes[n - 1] #max of the list
		poincaretimes = [] #Points for poincare section, x1 - x2 = pi
		#print np.array_str (taus)
		foutpop.write('{:.6f} {:.1f} {}'.format(T, lambda0timesTd, n))
		print np.array_str (photontimes)
		print 'n/T ' + str(n / T) + ' lambda0 ' + str(lambda0)

		t = 0 #Current time
		index = 0 #Awaiting the next photon

		#Random (ridiculous) conditions, so it can settle into normal range after a while.
		x1 = 10 * random.random()
		x2 = 10 * random.random()
		lastt = 0 #Last time since a photon was detected

		dec1 = np.exp(-1/T1*dt) #precompute the decline in one tick of x1,x2
		dec2 = np.exp(-1/T2*dt)
		probs = []

		timestart = time.clock()

		x1hist = [.7834] * int(Td / dt)
		x2hist = [0] * int(Td / dt)

		xdiff = 0

		while t < T:
			while index < n and photontimes[index] < t:
				#Get all the photons inside this tick
				ptime = photontimes[index]
				prob = np.sin( x1hist[0] - x2hist[0] + phi) ** 2

				if ptime > transcount * Td:
					probs.append(prob)
				#fout.write("{:.4e} {:.4e}\n".format(x1hist[0], x2hist[0]))

				#A photon is queued up, send it through the modulator
				if random.random() <= prob:
					#Photon passed through!
					#print 'Pop'
					if ptime > transcount * Td:
						#Filter out transients
						#photonpops.append(ptime)
						foutpop.write("{:7f}\n".format(ptime))
					x1 += beta / lambda0
					x2 += beta / lambda0

				index += 1

			x1 *= dec1
			x2 *= dec2
			#foutx.write("{:6f} {:6f}\n".format(x1, x2))

			#Take the poincare slice
			if (x1 - x2 - np.pi) * xdiff < 0:
				#poincaretimes.append(t)
				foutx.write("{:7f}\n".format(t - transcount * Td))
			xdiff = x1 - x2 - np.pi

			#Progress
			if index % (n / 50) == 0:
				print '=' * int(index / (n / 50)) + str(int(100.*index/n)) + "%"

			#Update history
			x1hist.append(x1)
			x2hist.append(x2)
			x1hist.pop(0)
			x2hist.pop(0)

			t += dt

		timeend = time.clock()

		print 'Simulation complete! Time elapsed = ' + str(timeend - timestart)

		foutpop.close()
		foutx.close()

		print 'Files saved!'
else:
	T = 5. #MUST BE FLOAT, lenght of sim time
	transtime = 0.1 #drop transients, so time to start counting data
	samptime = 0.00001
	sampspersec = 1 / samptime #number of samples to collect

	t = 0 #Current time
	#Random (ridiculous) conditions, so it can settle into normal range after a while.
	x1 = 10 * random.random()
	x2 = 10 * random.random()
	lastt = 0 #Last time since a photon was detected

	x1hist = [.7834] * int(Td / dt)
	x2hist = [0] * int(Td / dt)
	Ihist = [.123] * int(Td / dt)

	timestart = time.clock()

	fout = open("detint.out","w")
	foutv = open("detv.out","w")
	foutx = open("detxs.out","w")
	fout.write("{:6f}\n".format(T))
	fout.write("{}\n".format(sampspersec))

	intensities = []
	voltages = []
	timegraph = []
	poincaretimes = []
	xdiff = 0
	while t < T:
		I = (np.sin(x1hist[0] - x2hist[0] + phi)) ** 2
		x1 += (- 1 / T1 * x1 + beta * I) * dt
		x2 += (- 1 / T2 * x2 + beta * I) * dt

		if t > transtime and int(t / dt) % int(samptime / dt) == 0:
			#intensities.append(I)
			#voltages.append(x1 - x2)
			#timegraph.append(t - transtime)

			fout.write("{:6f}\n".format(I))
			foutv.write("{:6f}\n".format(x1 - x2))
			if (x1 - x2 - np.pi) * xdiff < 0: 
				foutx.write("{:7f}\n".format(t - transtime))
				poincaretimes.append(t - transtime)
			xdiff = (x1 - x2 - np.pi)

		#Progress
		if int(t / dt) % int(float(T) / 50 / dt) == 0:
			print int(t / dt) / int(float(T) / 100 / dt)

		#Update history
		x1hist.append(x1)
		x2hist.append(x2)
		x1hist.pop(0)
		x2hist.pop(0)
		t += dt
	fout.close()
	foutv.close()
	foutx.close()

	print 'Files saved!'

print 'Program done...'
