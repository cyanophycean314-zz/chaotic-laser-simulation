#Process experimental data to the vt.out format

import sys
import matplotlib.pyplot as plt
import numpy as np

filelist = [sys.argv[1]]

#The (basically unchangeable conditions)
dt = 0.00002 #Time interval (around the size of mu for lambda0 * Td = 3200)
samptime = 0.00002 #How often to take a sample
sampspersec = 1 / samptime #Inverse
Td = 0.001734 #Time delay
T1 = 0.0012 #Time constant for variable 1
T2 = 0.00006 #Time constant for variable 2
transtime = 0.1
phi = np.pi / 4 #Filter phase displacement
veps = 0.07
psliceval = 3.1415926535

#Simulation parameters
betatimesTd = 8.87 #this is the actual measurement that Aaron used, different than what he claims
beta = betatimesTd / Td #this is the real beta value, in the thousands.

for filename in filelist:
	fin = open(filename,'r')
	#lambda0 = float(raw_input("Photon rate (counts/s) = "))
	offset = float(sys.argv[2])#float(raw_input("Offset (V) = "))
	fout = open(filename[:-4] + 'vt.out','w')
	foutv = open(filename[:-4] + 'v.out','w')
	few = 3
	pastfew = [-1] * few
	counter = 0
	for i in range(3):
		fin.readline()
	t1 = float(fin.readline().split(',')[0])
	t2 = float(fin.readline().split(',')[0])
	deltat = t2 - t1
	N1 = max(1, int(Td / 4 / deltat))
	v1hist = [-1] * (2 * N1)
	vdiff = 0

	fout.write('-1\n')
	allts = []
	allvs = []
	lineno = 0
	for line in fin:
		 if len(line.split(',')) == 2 and lineno > 2:

		 	traw, vraw = line.split(',')
		 	v = np.pi * (float(vraw) - offset) / 2 / 1.55 #Scale it

		 	foutv.write("{:6f}\n".format(v))
		 	allvs.append(v)
		 	allts.append(float(traw))
		 	if (v - psliceval) * vdiff < 0 and vraw not in pastfew and len(allvs) > N1:
		 		fout.write("{:6f} {:6f}\n".format(v1hist[0], v1hist[N1]))
		 	vdiff = v - psliceval
		 	pastfew[counter] = vraw
		 	v1hist.append(v)
		 	v1hist.pop(0)
		 	counter = (counter + 1) % few
		 else:
		 	lineno += 1

	lengthofgraph = int(0.5/deltat)
	fout.close()
	plt.figure(1, figsize = (15,10))
	plt.title(filename)
	plt.subplot(211)
	plt.plot([x - allts[0] for x in allts[:lengthofgraph]], allvs[:lengthofgraph])	
	plt.subplot(212)
	plt.hist(allvs, bins = 0.6326 + 0.0924*np.arange(70))
	plt.savefig(filename[:-4] + ".png")
