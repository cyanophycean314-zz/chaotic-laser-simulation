# chaotic-laser-simulation

Exp_sim is the main program. Run simulations or read files. Easier to loop, has graph functions
Exp_sim_old is some old version of the program
hardvariancefinder is an old version of Exp_sim that tries to find the variance with the complicated theoretical formula
Kullbackleibler is for finding Kullback-leibler divergences
oldgraphy is an old graph program (functionalized ewww)

Data Batches
0 - First run.
	pop and x - Records photon pops and x1, x2
	
1 - Various photon rates from 12.5 to 3200, generally small. Also includes deterministic case
	pop and xs - Records photon pops and times for x1 - x2 = pi
	
2 - Various beta values from 0 to 9
	pop and xs - Records photon pops and times for x1 - x2 = pi
	
3 - Various photon rates from 1000 to 20000
	pop and xs - Records photon pops and times for x1 - x2 = pi
	
4 - Various photon rates from 10 to 30000
	pop and xs and v - Records photon pops and the times for x1 - x2 = pi and voltages
	
5 - Various photon rates from 1000 to 30000. also deterministic
	xs and v and vf - Records times for x1 - x2 = pi and voltages and smoothed voltages when x1 - x2 = pi (vf didn't work)
	
6 - Gaussian noise added to the deterministic case (model not the best)
	vt - Voltage at two time delays before, voltage at one time delay before
	
7 - 150000NT idk what
	vt - Voltage at two time delays before, voltage at one time delay before
	
8 - Various photon rates. Just archived cause lots of data to be added. Sucks we have different amounts of data
	vt - Voltage at two time delays before, voltage at one time delay before
	
9 - Experimental data before m.G and beta is really resolved.
	v,vt - Now examining time series again
