#Poincare generator

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
n = 100000000 #Photons to generate
xeps = 0.003 #X epsilon
eps = 0.00001 #other epsilon

simulation = "3200big"#"200" #Do the simulation or read a file
autocorr = False
lowlevel = False #Discrete or continuous simulation
looptime = True #Do time-steps or photon-steps

#Analysis begins
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

print 'Finished file reading'

w = Td/4
binw = w/100
binno = int(T / binw)

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

def movingwindow(t):
	#print binno
	lowerbound = binsearch(photonpops, t - w)
	upperbound = binsearch(photonpops, t)
	return upperbound - lowerbound

#Poincare section only necessary for Td/4
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
		nwp = movingwindow(ptime) #cause we already chopped off a transcount
		nwpt = movingwindow(ptime - Td / 4)
		if nwp >= maxp or nwpt >= maxp:
			#Don't wanna deal with outliers
			continue
		psec[int(nwp / pbinw)][int(nwpt / pbinw)] += 1

print 'Poincare section done!'

#3D Diagram
bigplot = [[],[],[]]
for i in range(100,20000):
	tt = transcount * Td + binw * i
	bigplot[0].append(movingwindow(tt))
	bigplot[1].append(movingwindow(tt - Td / 4))
	bigplot[2].append(movingwindow(tt - 2 * Td / 4))

print '3d attractor done!'

plt.figure(1)
plt.pcolor(np.array(psec))
fig = plt.figure(2)
ax = fig.add_subplot(111, projection='3d')
ax.plot(bigplot[0],bigplot[1],bigplot[2])
plt.show()