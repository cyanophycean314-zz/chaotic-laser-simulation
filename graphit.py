#Graph the files

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

#The (basically unchangeable conditions)
dt = 0.000005 #Time interval (around the size of mu for lambda0 * Td = 3200)
samptime = 0.0001 #How often to take a sample
sampspersec = 1 / samptime #Inverse
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
phi = np.pi / 4 #Filter phase displacement
#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.
T = 5 #seconds to simulate
deterministic = False

filelist = []
histogram = False
autocorr = False
poincare = True
attractor3d = False
variance = False

for filename in filelist:
	finv = open(str(filename) + "v.out","r")
	finx = open(str(filename) + "xs.out","r")

	poincaretimes = []
	for line in finx:
		poincaretimes.append(float(line))

	voltages = []
	for line in finv:
		voltages.append(float(line))
	T = voltages.pop(0)

	window = int(Td / 4 * sampspersec)
	smoothedvoltages = [np.mean(voltages[i - window: i]) for in range(window, len(voltages))]
	voltages = smoothedvoltages
	timegraph = [x * samptime for x in range(len(smoothedvoltages))]

	#Graph the stuff