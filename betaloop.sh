#!/bin/bash
#Run beta values

echo 'I run beta loop'

for i in 7 8 9; do
	for j in 0 1 2 3 4 5 6 7 8 9; do
		python exp_sim.py 5000 $i.$j
	done
done