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

simulation = "3200"

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

print 'Done reading!'

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
probno = 200

lowcnt = 0#binsearch(photonpops, transcount * Td * 0.99)
highcnt = 0#binsearch(photonpops, transcount * Td)
for i in range(offset, offset + binno):
	#Sliding count - linear, suggested by Joe
	while lowcnt < len(photonpops) and photonpops[lowcnt] < i * binw:
		lowcnt += 1
	while highcnt < len(photonpops) and photonpops[highcnt] < i * binw + w:
		highcnt += 1
	movingwindow[i - offset] = highcnt - lowcnt

mwhistprob = np.histogram(movingwindow, bins = probno)[0] / float(len(movingwindow))

#Compare with the uniform
######################################################
#Make the uniform
uniformprob = np.ones(len(mwhistprob)) / len(mwhistprob)

#Time to go integrate
kldiv = 0
for i in range(probno):
	kldiv += mwhistprob[i] * np.log(mwhistprob[i] / uniformprob[i])
print kldiv

plt.figure(1)
plt.subplot(211)
plt.plot(mwhistprob)
plt.subplot(212)
plt.plot(uniformprob)
plt.show()