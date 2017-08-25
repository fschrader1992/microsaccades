#!/bin/bash
#read video name

#here use poisson statistics 
for i in {1..10} 
do

    python ms_network.py poletti2010/exp1/cond3 exp1_cond3 121. 121. 
    python ms_network.py poletti2010/exp2/cond5 exp2_cond5 121. 121. 

done


#since there are already ten, no statistics from the network part
for i in {1..10} 
do

    python ms_network.py poletti2010/exp1/cond1/${i} exp1_cond1_${i} 121. 121. 
    python ms_network.py poletti2010/exp1/cond1_light/${i} exp1_cond1_light_${i} 121. 121. 
    python ms_network.py poletti2010/exp1/cond3_light/${i} exp1_cond3_light_${i} 121. 121. 
    python ms_network.py poletti2010/exp2/cond1/${i} exp2_cond1_${i} 121. 121. 
    python ms_network.py poletti2010/exp2/cond3/${i} exp2_cond3_${i} 121. 121. 

done

for a in {0..15} #15
do

    python ms_network.py poletti2010/exp1/cond2/${a}_8_pi_arc exp1_cond2_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp1/cond2_light/${a}_8_pi_arc exp1_cond2_light_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp1/cond4/${a}_8_pi_arc exp1_cond4_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc exp1_cond4_light_${a}_8_pi_arc 121. 121.
    
    python ms_network.py poletti2010/exp2/cond2/${a}_8_pi_arc exp2_cond2_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp2/cond4/${a}_8_pi_arc exp2_cond4_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp2/cond6/${a}_8_pi_arc exp2_cond6_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc exp2_cond7_${a}_8_pi_arc 121. 121.
    python ms_network.py poletti2010/exp2/cond8/${a}_8_pi_arc exp2_cond8_${a}_8_pi_arc 121. 121.

done