#!/bin/bash

date
echo 'STARTING SIMULATIONS'

for a in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_network_rest_weight.py poletti2010/exp1/cond2/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_${a}_8_pi_arc_${i} 280. 280. 1 2
        python ms_network_rest_weight.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_light_${a}_8_pi_arc_${i} 280. 280. 1 2_light
        python ms_network_rest_weight.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 280. 280. 1 4
        #python ms_network_rest_weight.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 
        python ms_network_rest_weight.py poletti2010/exp2/cond2/${a}_8_pi_arc/ ${i} ${i} exp2_cond2_${a}_8_pi_arc_${i} 280. 280. 2 2    
        python ms_network_rest_weight.py poletti2010/exp2/cond4/${a}_8_pi_arc/ ${i} ${i} exp2_cond4_${a}_8_pi_arc_${i} 280. 280. 2 4 
        python ms_network_rest_weight.py poletti2010/exp2/cond6/${a}_8_pi_arc/ ${i} ${i} exp2_cond6_${a}_8_pi_arc_${i} 280. 280. 2 6 
        #python ms_network_rest_weight.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 
        python ms_network_rest_weight.py poletti2010/exp2/cond8/${a}_8_pi_arc/ ${i} ${i} exp2_cond8_${a}_8_pi_arc_${i} 280. 280. 2 8 
        date
        echo ${a}: all files processed
    done
done

for i in {1..11} 
do

    date
    python ms_network_rest_weight.py poletti2010/exp1/cond1/ ${i} ${i} exp1_cond1_${i} 280. 280. 1 1
    python ms_network_rest_weight.py poletti2010/exp1/cond1_light/ ${i} ${i} exp1_cond1_light_${i} 280. 280. 1 1_light
    python ms_network_rest_weight.py poletti2010/exp1/cond3/ 1 ${i} exp1_cond3 280. 280. 1 3
    python ms_network_rest_weight.py poletti2010/exp1/cond3_light/ ${i} ${i} exp1_cond3_light_${i} 280. 280. 1 3_light 
    python ms_network_rest_weight.py poletti2010/exp2/cond1/ ${i} ${i} exp2_cond1_${i} 280. 280. 2 1
    python ms_network_rest_weight.py poletti2010/exp2/cond3/ ${i} ${i} exp2_cond3_${i} 280. 280. 2 3
    python ms_network_rest_weight.py poletti2010/exp2/cond5/ 1 ${i} exp2_cond5 280. 280. 2 5
    date
    echo ${i}: five files processed
done

echo 'exp 1 con 1, 1_light, 3 and 3 simulated'

for a in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 #1 2 3 5 7 10 #0 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date

        python ms_network_rest_weight.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 280. 280. 1 4_light
        python ms_network_rest_weight.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 280. 280. 2 7
        date
        echo ${a}: all files processed
    done
done

date
echo done with this


#!/bin/bash

date
echo 'STARTING SIMULATIONS'

for a in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_network_rest_modet_weight.py poletti2010/exp1/cond2/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_${a}_8_pi_arc_${i} 280. 280. 1 2
        python ms_network_rest_modet_weight.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_light_${a}_8_pi_arc_${i} 280. 280. 1 2_light
        python ms_network_rest_modet_weight.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 280. 280. 1 4
        #python ms_network_rest_modet_weight.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 
        python ms_network_rest_modet_weight.py poletti2010/exp2/cond2/${a}_8_pi_arc/ ${i} ${i} exp2_cond2_${a}_8_pi_arc_${i} 280. 280. 2 2    
        python ms_network_rest_modet_weight.py poletti2010/exp2/cond4/${a}_8_pi_arc/ ${i} ${i} exp2_cond4_${a}_8_pi_arc_${i} 280. 280. 2 4 
        python ms_network_rest_modet_weight.py poletti2010/exp2/cond6/${a}_8_pi_arc/ ${i} ${i} exp2_cond6_${a}_8_pi_arc_${i} 280. 280. 2 6 
        #python ms_network_rest_modet_weight.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 
        python ms_network_rest_modet_weight.py poletti2010/exp2/cond8/${a}_8_pi_arc/ ${i} ${i} exp2_cond8_${a}_8_pi_arc_${i} 280. 280. 2 8 
        date
        echo ${a}: all files processed
    done
done

for i in {1..11} 
do

    date
    python ms_network_rest_modet_weight.py poletti2010/exp1/cond1/ ${i} ${i} exp1_cond1_${i} 280. 280. 1 1
    python ms_network_rest_modet_weight.py poletti2010/exp1/cond1_light/ ${i} ${i} exp1_cond1_light_${i} 280. 280. 1 1_light
    python ms_network_rest_modet_weight.py poletti2010/exp1/cond3/ 1 ${i} exp1_cond3 280. 280. 1 3
    python ms_network_rest_modet_weight.py poletti2010/exp1/cond3_light/ ${i} ${i} exp1_cond3_light_${i} 280. 280. 1 3_light 
    python ms_network_rest_modet_weight.py poletti2010/exp2/cond1/ ${i} ${i} exp2_cond1_${i} 280. 280. 2 1
    python ms_network_rest_modet_weight.py poletti2010/exp2/cond3/ ${i} ${i} exp2_cond3_${i} 280. 280. 2 3
    python ms_network_rest_modet_weight.py poletti2010/exp2/cond5/ 1 ${i} exp2_cond5 280. 280. 2 5
    date
    echo ${i}: five files processed
done

echo 'exp 1 con 1, 1_light, 3 and 3 simulated'

for a in  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 #1 2 3 5 7 10 #0 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date

        python ms_network_rest_modet_weight.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 280. 280. 1 4_light
        python ms_network_rest_modet_weight.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 280. 280. 2 7
        date
        echo ${a}: all files processed
    done
done

date
echo done with this

