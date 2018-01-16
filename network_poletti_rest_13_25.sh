#!/bin/bash
for i in {1..10} 
do

    date
    python ms_network_rest.py poletti2010/exp1/cond3/ 1 ${i} exp1_cond3 280. 280. 1 3
    python ms_network_rest.py poletti2010/exp2/cond5/ 1 ${i} exp2_cond5 280. 280. 2 5
    date
    echo ${i}: five files processed
done