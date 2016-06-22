#!/usr/bin/env python
# Simulate the differential equations in Hagerstrom (2015)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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

#Conditions
dt = 0.0000005 #Time interval (around the size of mu for lambda0 * Td = 3200)
Td = 0.001734 #Time delay
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
lambda0timesTd = 600 #Metric given in Aaron's paper
lambda0 = lambda0timesTd / Td
mu = 1 / lambda0 #Poisson interarrival time average
phi = np.pi / 4 #Filter phase displacement
transcount = 100 #How many Td's to wait for transients to die
n = 6000000 #Photons to generate
xeps = 0.003 #X epsilon
eps = 0.00001 #other epsilon

simulation = "3200"#"200" #Do the simulation or read a file
autocorr = False #Do the autocorrelation, kinda slow
poincare = False
variance = False

if simulation == "simulation":
	foutpop = open('pop.out','w')
	foutx = open('xs.out','w') #x-simplified

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
else:
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

print 'Now graphing...'

#Print results
######################################


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

w = Td / 4
binw = w / 20
offset = int(transcount * Td / binw)
binno = int(T / binw) - offset
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

#Autocorrelation
#############################
if autocorr:
	tinterval = 2 * offset #length of the section of the graph we will be shifting and comparing
	mws = movingwindow[:tinterval] #mws = movingwindowslice
	iterlen = offset / 20 #how far to shift along the graph
	mean = np.average(mws)
	C = []

	for shift in range (iterlen):
		j = shift
		ans = 0
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
else:
	Cgraph = [1]

#Poincare section
################################
if poincare:
	if lambda0timesTd < 13:
		pbinw = 1
		maxp = 8
	elif lambda0timesTd < 210:
		pbinw = 1
		maxp = 80
	else:
		pbinw = 4
		maxp = 800

	psec = [[int(_) for _ in x] for x in np.zeros((maxp / pbinw, maxp / pbinw))] #(N_w(t), N_w(t - Td/4))
	for ptime in poincaretimes:
		if ptime >= transcount * Td + Td / 4:
			nwp = movingwindowcalc(ptime) #cause we already chopped off a transcount
			nwpt = movingwindowcalc(ptime - Td / 4)
			if nwp >= maxp or nwpt >= maxp:
				#Don't wanna deal with outliers
				continue
			psec[int(nwp / pbinw)][int(nwpt / pbinw)] += 1

	print 'Poincare section done!'

	#3D Diagram
	bigplot = [[],[],[]]
	for i in range(100,20000):
		tt = transcount * Td + binw * i
		bigplot[0].append(movingwindowcalc(tt))
		bigplot[1].append(movingwindowcalc(tt - Td / 4))
		bigplot[2].append(movingwindowcalc(tt - 2 * Td / 4))

	print '3d attractor done!'

#Variance
#############################
if variance:
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
		for i in range(offset, offset + binno):
			lowerbound = binsearch(photonpops, binw * i)
			upperbound = binsearch(photonpops, binw * i + w) #Don't start with zero counts
			mws[i - offset] = upperbound - lowerbound #no +1 because binsearch returns number over the key

		varOverW.append(np.var(mws) / w)

#Poincare Slice
####################################

plt.figure(1)
plt.title('lambda0 * Td = ' + str(lambda0timesTd))
plt.subplot(311)
plt.xlim([0,T - transcount * Td])
plt.plot(timegraph, movingwindow)
plt.subplot(312)
plt.hist(movingwindow, bins = range(0,900,3))
plt.subplot(313)
plt.plot(Cgraph)

if poincare:
	plt.figure(2)
	plt.pcolor(np.array(psec))
	fig = plt.figure(3)
	ax = fig.add_subplot(111, projection='3d')
	ax.plot(bigplot[0],bigplot[1],bigplot[2])

if variance:
	plt.figure(4)
	plt.loglog(Ws, varOverW)

plt.show()