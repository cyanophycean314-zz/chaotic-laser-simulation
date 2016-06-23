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
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
phi = np.pi / 4 #Filter phase displacement
xeps = 0.003 #X epsilon
eps = 0.00001 #other epsilon
#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
lambda0timesTd = 3200 #Metric given in Aaron's paper
lambda0 = lambda0timesTd / Td
mu = 1 / lambda0 #Poisson interarrival time average
transcount = 100 #How many Td's to wait for transients to die
n = 5000000 #Photons to generate
#Histogram parameters
w = Td / 4
binw = w / 20
offset = int(transcount * Td / binw)

simulation = True
deterministic = False
#filelist = ["lambda0Td=" + str(x) for x in [1000,1500,2000,2500,3000,3500,4000,5000,10000,20000]]
#filelist = ["lambda0Td=10000"]
histogram = True
autocorr = False
poincare = True
attractor3d = False
variance = False

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
	lowerbound = binsearch(photonpops, t - w)
	upperbound = binsearch(photonpops, t)
	return upperbound - lowerbound

#Generate graphs
#Get histograms
####################################
def getmovingwindow(photonpops):
	binno = int(T / binw) - offset + 1
	movingwindow = [0] * (binno)
	#print binno

	lowcnt = 0 #binsearch(photonpops, transcount * Td * 0.99)
	highcnt = 0 #binsearch(photonpops, transcount * Td)
	for i in range(offset, offset + binno):
		#Sliding count - linear, suggested by Joe
		while lowcnt < len(photonpops) and photonpops[lowcnt] < i * binw:
			lowcnt += 1
		while highcnt < len(photonpops) and photonpops[highcnt] < i * binw + w:
			highcnt += 1
		movingwindow[i - offset] = highcnt - lowcnt

	#print movingwindow
	#Toss out initial transient phase
	print 'Histogram complete!'
	timegraph = [x * binw for x in range(binno)]

	return timegraph, movingwindow

#Autocorrelation
#############################
def getautocorr(Nw, deterministic=False):
	if not deterministic:
		tinterval = 2 * offset #length of the section of the graph we will be shifting and comparing
		mws = Nw[:tinterval] #mws = movingwindowslice
		iterlen = offset / 20 #how far to shift along the graph
	else:
		tinterval = int(0.05 * sampspersec)
		mws = Nw[:tinterval]
		iterlen = int(0.02 * sampspersec)

	mean = np.average(mws)
	C = []
	for shift in range (iterlen):
		j = shift
		ans = 0.
		while j < tinterval:
			ans += (mws[j-shift] - mean) * (mws[j] - mean)
			j += 1
		ans /= (tinterval - shift)
		C.append(ans)

		#Progress
		if shift % (iterlen / 10) == 0:
			print shift * 100. / iterlen

	C = [x / C[0] for x in C] #Normalize it

	#Create the symmetric C
	Cgraph = C[::-1] + [1] + C
	#print Cgraph
	print 'Autocorrelation done!'
	return Cgraph

#Poincare section
################################
def getpoincare(Nw, poincaretimes, ltTd, deterministic=False):
	if not deterministic:
		if ltTd < 13:
			pbinw = 1
			maxp = 8
		elif ltTd < 210:
			pbinw = 1
			maxp = 80
		elif ltTd < 5000:
			pbinw = 4
			maxp = int(math.ceil(ltTd / 4 / pbinw) * pbinw)
		else:
			pbinw = 10
			maxp = int(math.ceil(ltTd / 4 / pbinw) * pbinw)

		movingwindow = Nw
		psec = [[int(_) for _ in x] for x in np.zeros((maxp / pbinw, maxp / pbinw))] #(N_w(t), N_w(t - Td/4))
		pslice = [0 for x in range(maxp / 2 / pbinw)]
		#The slice we're taking is y = maxp -2(x - maxp / 4)
		for ptime in poincaretimes:
			if ptime >= transcount * Td + Td / 4:
				#print str(int(ptime / binw) - offset) + "," + str(len(movingwindow))
				nwp = movingwindowcalc(ptime) #cause we already chopped off a transcount
				nwpt = movingwindowcalc(ptime - Td / 4)
				if nwp >= maxp or nwpt >= maxp:
					#Don't wanna deal with outliers
					continue
				#print str(nwp) + "," + str(nwpt) + "," + str(maxp)
				psec[int(nwp / pbinw)][int(nwpt / pbinw)] += 1
				if int((maxp - 2*(nwp - maxp / 4)) / pbinw) == int(nwpt / pbinw):
					print str(nwp) + "," + str(nwpt)
					pslice[int(nwp - maxp / 4) / pbinw] += 1

		print 'Poincare section and slice done!'
	else:
		intensities = Nw
		pbinw = 0.005
		maxp = 1
		psec = [[int(_) for _ in x] for x in np.zeros((maxp / pbinw, maxp / pbinw))] #(N_w(t), N_w(t - Td/4))
		pslice = [0 for x in range(int(maxp / pbinw))]
		for ptime in poincaretimes:
			if ptime > Td:
				nwp = intensities[int((ptime - transtime) * sampspersec)]
				nwpt = intensities[int((ptime - transtime - Td / 4) * sampspersec)]
				psec[int(nwp / pbinw)][int(nwpt / pbinw)] += 1
				pslice[int(nwp / pbinw)] += 1

		print 'Poincare section and slice done!'

	return psec, pslice, range(len(pslice))

def get3dattractor(Nw, deterministic=False):
	#3D Diagram
	bigplot = [[],[],[]]
	if deterministic:
		for i in range(int(Td * sampspersec / 2), 20000):
			tt = transcount * Td + binw * i
			bigplot[0].append(movingwindowcalc(tt))
			bigplot[1].append(movingwindowcalc(tt - Td / 4))
			bigplot[2].append(movingwindowcalc(tt - 2 * Td / 4))
	else:
		for i in range(int(2 * Td / 4 / binw),20000):
			bigplot[0].append(Nw[i])
			bigplot[1].append(Nw[i - int(Td / 4 / binw)])
			bigplot[2].append(Nw[i - int(2 * Td / 4 / binw)])

	print '3d attractor done!'
	return bigplot
#Variance
#############################
def getvariance(photonpops):
	varOverW = []
	Ws = []
	for wvary in list(range(35)):
		w = 10**(wvary / 5. - 5) * Td #Nice range for w
		binw = w / 40 #So we can calculate variance
		binno = transcount / 2 * int (w / binw) + max(1, transcount *int(2 * -np.log(binw)))
		Ws.append(w)
		offset = offset
		mws = [0] * (binno)
		counter = 0
		lowcnt = 0 #binsearch(photonpops, transcount * Td * 0.99)
		highcnt = 0 #binsearch(photonpops, transcount * Td)
		for i in range(offset, offset + binno):
			#Sliding count - linear, suggested by Joe
			while lowcnt < len(photonpops) and photonpops[lowcnt] < i * binw:
				lowcnt += 1
			while highcnt < len(photonpops) and photonpops[highcnt] < i * binw + w:
				highcnt += 1
			movingwindow[i - offset] = highcnt - lowcnt

		varOverW.append(np.var(mws) / w)

	return Ws, varOverW

#Poincare Slice
####################################
def showgraphs(Nw, timegraph, poincaretimes, ltTd, deterministic):
	if histogram:
		plt.figure(1)
		if deterministic:
			plt.subplot(311)
			plt.title('Deterministic')
			plt.xlim([0,T - transtime])
			plt.plot(timegraph, intensities)
			plt.subplot(312)
			plt.hist(intensities, bins = 20)
		else:
			plt.subplot(311)
			plt.title('lambda0 * Td = ' + str(lambda0timesTd))
			plt.xlim([0,T - transcount * Td])
			plt.plot(timegraph, Nw)
			plt.subplot(312)
			plt.hist(movingwindow, bins = 30)

		if autocorr:
			plt.subplot(313)
			Cgraph = getautocorr(Nw, deterministic)
			plt.plot(Cgraph)

	if poincare:
		psec, pslice, pbounds = getpoincare(Nw, poincaretimes, ltTd, deterministic)
		plt.figure(2)
		plt.subplot(211)
		plt.title("lambda0Td = {}".format(int(ltTd)))
		plt.pcolor(np.array(psec))
		plt.plot([maxp / 4, maxp * 3 / 4], [maxp, 0], color ='r')
		plt.subplot(212)
		plt.scatter(pbounds, pslice)

	if attractor3d:
		bigplot = get3dattractor(Nw, deterministic)
		fig = plt.figure(3)
		ax = fig.add_subplot(111, projection='3d')
		ax.plot(bigplot[0],bigplot[1],bigplot[2])

	if variance:
		Ws, varOverW = getvariance(photonpops)
		plt.figure(4)
		plt.loglog(Ws, varOverW)

	plt.show()

if simulation:
	if not deterministic:
		for lambda0timesTd in [1000,1500,2000,2500,3000,3500,4000,5000,10000,20000]:
			print lambda0timesTd
			lambda0 = lambda0timesTd / Td
			mu = 1 / lambda0 #Poisson interarrival time average
			n = 10000 * lambda0timesTd

			foutpop = open('lambda0Td=' + str(lambda0timesTd) + 'pop.out','w')
			foutx = open('lambda0Td=' + str(lambda0timesTd) + 'xs.out','w') #x-simplified
			#Generate photon times!
			taus = np.random.exponential (mu , n) #Interarrival times
			photontimes = np.cumsum(taus) #Cumulative sum, times of arrival at modulator
			photonpops = [] #Times of photon arivals at counter!
			T = photontimes[n - 1] #max of the list
			poincaretimes = [] #Points for poincare section, x1 - x2 = pi
			#print np.array_str (taus)
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
							photonpops.append(ptime)
						x1 += beta / lambda0
						x2 += beta / lambda0

					index += 1

				x1 *= dec1
				x2 *= dec2
				#foutx.write("{:6f} {:6f}\n".format(x1, x2))

				if abs(x1 - x2 - np.pi) < xeps:
					poincaretimes.append(t)
					foutx.write("{:7f}\n".format(t))

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

			#Output and save results
			foutpop.write('{:.6f} {:.1f} {}'.format(T, lambda0timesTd, n))
			for i in range(len(photonpops)):
				foutpop.write(str('{:.7f}\n'.format(photonpops[i])))
				if i % (len(photonpops) / 20) == 0:
					print i * 100 / len(photonpops) + 1

			foutpop.close()
			foutx.close()

			print 'Files saved!'
	else:
		T = 5 #MUST BE FLOAT, lenght of sim time
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

		timestart = time.clock()

		fout = open("detint.out","w")
		foutx = open("detxs.out","w")
		fout.write("{:6f}\n".format(T))
		fout.write("{}\n".format(sampspersec))

		intensities = []
		timegraph = []
		poincaretimes = []
		while t < T:
			x1 += (- 1 / T1 * x1 + beta * np.sin(x1hist[0] - x2hist[0] + phi) ** 2) * dt
			x2 += (- 1 / T2 * x2 + beta * np.sin(x1hist[0] - x2hist[0] + phi) ** 2) * dt

			if t > transtime and int(t / dt) % int(samptime / dt) == 0:
				intensities.append(np.sin(x1hist[0] - x2hist[0] + phi) ** 2)
				timegraph.append(t - transtime)
				fout.write("{:6f}\n".format(np.sin(x1hist[0] - x2hist[0] + phi) ** 2))
				if abs(x1 - x2 - np.pi) < xeps:
					foutx.write("{:7f}\n".format(t))
					poincaretimes.append(t)

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

		print 'Files saved!'
else:
	if not deterministic:
		for filename in filelist:
			finpop = open(str(filename) + "pop.out", 'r')
			finxs = open(str(filename) + "xs.out", 'r')
			photonpops = []
			poincaretimes = []
			for line in finpop:
				photonpops.append(line)
			photonpops.pop() #Last line is a bit buggy usually
			runinfo = photonpops.pop(0).split() #First line is the info sectioN
			for line in finxs:
				poincaretimes.append(float(line))
			for i in range(len(photonpops)):
				photonpops[i] = float(photonpops[i])
			n = int(float(runinfo[2]))
			lambda0timesTd = int(float(runinfo[1]))
			T = float(runinfo[0])

			#Graph
			print filename
			if histogram:
				timegraph, movingwindow = getmovingwindow(photonpops)
			else:
				timegraph = []
				movingwindow = []
			showgraphs(movingwindow, timegraph, poincaretimes, lambda0timesTd, deterministic)
	else:
		transtime = 0.1
		fin = open('detint.out','r')
		finxs = open('detxs.out','r')
		intensities = []
		for line in fin:
			intensities.append(float(line))
		T = intensities.pop(0)
		sampspersec = intensities.pop(0)
		timegraph = [x / sampspersec for x in range(int((T - transtime) * sampspersec))]

		poincaretimes = []
		for line in finxs:
			poincaretimes.append(float(line))
print 'Program done...'