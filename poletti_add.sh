#!/bin/bash

#python ms_input.py poletti2010/exp1/cond1/11 exp1_cond1_11 1000 &
#python ms_input.py poletti2010/exp2/cond1/11 exp2_cond1_11 1000 &
#python ms_input.py poletti2010/exp1/cond1/8 exp1_cond1_8 1000 &
#python ms_input.py poletti2010/exp2/cond1/8 exp2_cond1_8 1000 
#wait

python ms_input.py poletti2010/exp1/cond1/9 exp1_cond1_9 1000 &
python ms_input.py poletti2010/exp2/cond1/9 exp2_cond1_9 1000 &
python ms_input.py poletti2010/exp1/cond1/10 exp1_cond1_10 1000 &
python ms_input.py poletti2010/exp2/cond1/10 exp2_cond1_10 1000 
wait

#python ms_input.py poletti2010/exp1/cond4/1_8_pi_arc/1 exp1_cond4_1_8_pi_arc_1 1000 &
#python ms_input.py poletti2010/exp1/cond4/5_8_pi_arc/1 exp1_cond4_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond1/11 exp1_cond1_11 1000 &
python ms_input.py poletti2010/exp1/cond4/7_8_pi_arc/1 exp1_cond4_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4/12_8_pi_arc/1 exp1_cond4_12_8_pi_arc_1 1000 
wait

date
echo done with input

for a in 0 2 8 13 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp1/cond2/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond2/${a}_8_pi_arc/network/${i};
        fi
        #if [ ! -d ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc/network/${i} ]; then
        #    mkdir -p ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc/network/${i};
        #fi
        if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond2/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond2/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond4/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond4/${a}_8_pi_arc/network/${i};
        fi
        #if [ ! -d ./data/poletti2010/exp2/cond6/${a}_8_pi_arc/network/${i} ]; then
        #    mkdir -p ./data/poletti2010/exp2/cond6/${a}_8_pi_arc/network/${i};
        #fi
        #if [ ! -d ./data/poletti2010/exp2/cond8/${a}_8_pi_arc/network/${i} ]; then
        #    mkdir -p ./data/poletti2010/exp2/cond8/${a}_8_pi_arc/network/${i};
        #fi
        python ms_network.py poletti2010/exp1/cond2/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_${a}_8_pi_arc_${i} 241. 241. 1 2
        #python ms_network.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_light_${a}_8_pi_arc_${i} 241. 241. 1 2_light
        python ms_network.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 241. 241. 1 4
        #python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 
        python ms_network.py poletti2010/exp2/cond2/${a}_8_pi_arc/ ${i} ${i} exp2_cond2_${a}_8_pi_arc_${i} 241. 241. 2 2    
        python ms_network.py poletti2010/exp2/cond4/${a}_8_pi_arc/ ${i} ${i} exp2_cond4_${a}_8_pi_arc_${i} 241. 241. 2 4 
        #python ms_network.py poletti2010/exp2/cond6/${a}_8_pi_arc/ ${i} ${i} exp2_cond6_${a}_8_pi_arc_${i} 241. 241. 2 6 
        #python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 
        #python ms_network.py poletti2010/exp2/cond8/${a}_8_pi_arc/ ${i} ${i} exp2_cond8_${a}_8_pi_arc_${i} 281. 281. 2 8 
        date
        echo ${a}: all files processed
    done
done

for a in 1 5 7 12 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i};
        fi
        python ms_network.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 241. 241. 1 4
         date
        echo ${a}: all files processed
    done
done

for i in 8 9 10 11 
do
    if [ ! -d ./data/poletti2010/exp1/cond1/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp2/cond1/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond1/network/${i};
    fi
    date
    python ms_network.py poletti2010/exp1/cond1/ ${i} ${i} exp1_cond1_${i} 241. 241. 1 1
    python ms_network.py poletti2010/exp2/cond1/ ${i} ${i} exp2_cond1_${i} 241. 241. 2 1
    date
    echo ${i}: five files processed
done

echo network rest start


for a in 0 2 8 13 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_network_rest.py poletti2010/exp1/cond2/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_${a}_8_pi_arc_${i} 241. 241. 1 2
        #python ms_network.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_light_${a}_8_pi_arc_${i} 241. 241. 1 2_light
        python ms_network_rest.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 241. 241. 1 4
        #python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 
        python ms_network_rest.py poletti2010/exp2/cond2/${a}_8_pi_arc/ ${i} ${i} exp2_cond2_${a}_8_pi_arc_${i} 241. 241. 2 2    
        python ms_network_rest.py poletti2010/exp2/cond4/${a}_8_pi_arc/ ${i} ${i} exp2_cond4_${a}_8_pi_arc_${i} 241. 241. 2 4 
        #python ms_network.py poletti2010/exp2/cond6/${a}_8_pi_arc/ ${i} ${i} exp2_cond6_${a}_8_pi_arc_${i} 241. 241. 2 6 
        #python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 
        #python ms_network.py poletti2010/exp2/cond8/${a}_8_pi_arc/ ${i} ${i} exp2_cond8_${a}_8_pi_arc_${i} 281. 281. 2 8 
        date
        echo ${a}: all files processed
    done
done

for a in 1 5 7 12 #1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_network_rest.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 241. 241. 1 4
        date
        echo ${a}: all files processed
    done
done

for i in 8 9 10 11 
do
    date
    python ms_network_rest.py poletti2010/exp1/cond1/ ${i} ${i} exp1_cond1_${i} 241. 241. 1 1
    python ms_network_rest.py poletti2010/exp2/cond1/ ${i} ${i} exp2_cond1_${i} 241. 241. 2 1
    date
    echo ${i}: five files processed
done
