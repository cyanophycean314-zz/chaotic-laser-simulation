#!/bin/bash
#Run beta values

echo 'I run beta loop'

for i in 9; do
	for j in 0 1 2 3 4 5 6 7 8 9 10; do
		python exp_sim_noisytom.py 5000 $j
	done
done