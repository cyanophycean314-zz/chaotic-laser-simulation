#!/usr/bin/env python
# Simulate the differential equations in Hagerstrom (2015)

import numpy as np
import random
import time


'''
Outline

Generate photon times
Loop
	Get a photon
	Go through modulator
	Go through filter
'''

#The (basically unchangeable conditions)
dt = 0.00002 #Time interval (around the size of mu for lambda0 * Td = 3200)
samptime = 0.00002 #How often to take a sample
sampspersec = 1 / samptime #Inverse
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
transtime = 0.1
phi = np.pi / 4 #Filter phase displacement
xeps = 0.001

#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
deterministic = True
points = False #Count pops
T = 200. #seconds to simulate

if not deterministic:
	filelist = [30000]
	subscripts = [''] + list("abcdefghijklmno")#Allows for multiple files of the same photon rate
else:
	filelist = ["detx"]
	subscripts = ["005","008","01","03","05","07","1"]
	noises = [0.005, 0.008, 0.01, 0.03, 0.05, 0.07, 0.1]

for filename in filelist:
	for lettno in range(len(subscripts)):
		lett = subscripts[lettno]
		noise = noises[lettno]
		t = 0
		x1 = 10 * random.random()
		x2 = 10 * random.random()
		N = int(Td / dt)
		x1hist = [.7834] * N
		x2hist = [.7834] * N
		xdiff = 0
		pval = np.pi
		ctr = 0
		ctr2 = 0
		N2 = int(Td / samptime / 4)
		vhist = [0] * N2 #raw voltage over last Td / 4
		vhistsum = 0
		vlag = [2] * N2 #what the smoothed voltage was one Td / 4 ago
		N3 = 2 * N2
		vhisttwo = [0] * N3

		foutvt = open(str(filename) + lett + "vt.out","w")
		print str(filename) + lett
		#foutx = open(str(filename) + lett + "xs.out","w")

		timestart = time.clock()

		if not deterministic:
			if points:
				foutpop = open(str(filename) + "pop.out","w")
			lambda0timesTd = int(filename) #Metric given in Aaron's paper
			lambda0 = lambda0timesTd / Td
			mu = 1 / lambda0 #Poisson interarrival time average
			chunks = 100
			chunklen = int(T / mu / chunks)
			n = chunks * chunklen
			chunkno = 0
			taus = np.random.exponential(mu, chunklen)
			index = 0 #photon index
			lastt = 0
			dec1 = np.exp(-dt/T1)
			dec2 = np.exp(-dt/T2)

		foutvt.write(str(T) + "\n") #only predicted, you never know for sure :O

		while (deterministic and t < T) or (not deterministic and chunkno < chunks):
			I = (np.sin(x1hist[ctr % N] - x2hist[ctr % N] + phi)) ** 2

			#Evolution of x1, x2
			if deterministic:
				x1 += (-1 / T1 * x1 + beta * I) * dt
				x2 += (-1 / T2 * x2 + beta * I) * dt
			else:
				while index < len(taus) and t > lastt + taus[index]:
					if random.random() <= I:
						x1 += beta / lambda0
						x2 += beta / lambda0
						if points and t >= transtime:
							foutpop.write("{:6f}\n".format(t))
					lastt += taus[index]
					index += 1
					if index >= len(taus):
						index = 0
						taus = np.random.exponential(mu, chunklen)
						chunkno += 1

				x1 *= dec1
				x2 *= dec2

			#Record data
			if int(t / dt) % int(samptime / dt) == 0:
				v = x1 - x2
				if deterministic:
					v *= np.random.normal(1, noise) #Add the noise
				if t >= transtime:
					#foutv.write("{:6f}\n".format(x1 - x2))
					if (v - pval) * xdiff < 0:
						#foutx.write("{:6f}\n".format(t))
						foutvt.write("{:6f} {:6f}\n".format(vhisttwo[0], vhisttwo[N2]))
				vhisttwo.append(v)
				vhisttwo.pop(0)
				xdiff = v - pval
				ctr2 += 1

			#Progress
			if int(t / dt) % int(T / 50. / dt) == 0:
				percent = int(t / dt) / int(T / 100. / dt)
				print ('=' * (percent / 2)) + str(percent)

			x1hist[ctr % N] = x1
			x2hist[ctr % N] = x2
			t += dt
			ctr += 1

		#foutv.close()
		#foutx.close()
		foutvt.close()
		#foutvf.close()
		if points:
			foutpop.close()

		print str(filename) + ": " + str(time.clock() - timestart)

print 'Program done...'