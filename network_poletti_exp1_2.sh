#!/bin/bash

date
echo 'STARTING SIMULATIONS'

for a in 1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp1/cond2/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond2/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond2/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond2/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond4/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond4/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond6/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond6/${a}_8_pi_arc/network/${i};
        fi
        if [ ! -d ./data/poletti2010/exp2/cond8/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond8/${a}_8_pi_arc/network/${i};
        fi
        python ms_network.py poletti2010/exp1/cond2/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_${a}_8_pi_arc_${i} 241. 241. 1 2
        python ms_network.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond2_light_${a}_8_pi_arc_${i} 241. 241. 1 2_light
        python ms_network.py poletti2010/exp1/cond4/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_${a}_8_pi_arc_${i} 241. 241. 1 4
        #python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 
        python ms_network.py poletti2010/exp2/cond2/${a}_8_pi_arc/ ${i} ${i} exp2_cond2_${a}_8_pi_arc_${i} 241. 241. 2 2    
        python ms_network.py poletti2010/exp2/cond4/${a}_8_pi_arc/ ${i} ${i} exp2_cond4_${a}_8_pi_arc_${i} 241. 241. 2 4 
        python ms_network.py poletti2010/exp2/cond6/${a}_8_pi_arc/ ${i} ${i} exp2_cond6_${a}_8_pi_arc_${i} 241. 241. 2 6 
        #python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 
        python ms_network.py poletti2010/exp2/cond8/${a}_8_pi_arc/ ${i} ${i} exp2_cond8_${a}_8_pi_arc_${i} 281. 281. 2 8 
        date
        echo ${a}: all files processed
    done
done

for i in {2..7} 
do
    if [ ! -d ./data/poletti2010/exp1/cond1/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1/network/${i};
    fi
    
    if [ ! -d ./data/poletti2010/exp1/cond1_light/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1_light/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp1/cond3/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond3/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp1/cond3_light/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond3_light/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp2/cond1/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond1/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp2/cond3/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond3/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp2/cond5/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond5/network/${i};
    fi

    date
    python ms_network.py poletti2010/exp1/cond1/ ${i} ${i} exp1_cond1_${i} 241. 241. 1 1
    python ms_network.py poletti2010/exp1/cond1_light/ ${i} ${i} exp1_cond1_light_${i} 241. 241. 1 1_light
    python ms_network.py poletti2010/exp1/cond3/ 1 ${i} exp1_cond3_${i} 241. 241. 1 3
    python ms_network.py poletti2010/exp1/cond3_light/ ${i} ${i} exp1_cond3_light_${i} 241. 241. 1 3_light 
    python ms_network.py poletti2010/exp2/cond1/ ${i} ${i} exp2_cond1_${i} 241. 241. 2 1
    python ms_network.py poletti2010/exp2/cond3/ ${i} ${i} exp2_cond3_${i} 241. 241. 2 3
    python ms_network.py poletti2010/exp2/cond5/ 1 ${i} exp2_cond5_${i} 241. 241. 2 5
    date
    echo ${i}: five files processed
done

echo 'exp 1 con 1, 1_light, 3 and 3 simulated'

for a in  1 2 3 #0 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc/network/${i};
        fi

        if [ ! -d ./data/poletti2010/exp2/cond7/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond7/${a}_8_pi_arc/network/${i};
        fi

        python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 241. 241. 1 4_light
        python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 281. 281. 2 7
        date
        echo ${a}: all files processed
    done
done

date
echo done with this


python ms_input.py poletti2010/exp2/cond7/5_8_pi_arc/1 exp2_cond7_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/7_8_pi_arc/1 exp2_cond7_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/10_8_pi_arc/1 exp2_cond7_10_8_pi_arc_1 1000 
wait
date
echo third three files processed

for a in 5 7 10 #12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp2/cond7/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp2/cond7/${a}_8_pi_arc/network/${i};
        fi
        python ms_network.py poletti2010/exp2/cond7/${a}_8_pi_arc/ ${i} ${i} exp2_cond7_${a}_8_pi_arc_${i} 281. 281. 2 7
        date
        echo ${a}: all files processed
    done
done

python ms_input.py poletti2010/exp1/cond4_light/5_8_pi_arc/1 exp1_cond4_light_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/7_8_pi_arc/1 exp1_cond4_light_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/10_8_pi_arc/1 exp1_cond4_light_10_8_pi_arc_1 1000 
wait
date
echo fourth three files processed

for a in   5 7 10 #12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        if [ ! -d ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc/network/${i} ]; then
            mkdir -p ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc/network/${i};
        fi

        python ms_network.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/ ${i} ${i} exp1_cond4_light_${a}_8_pi_arc_${i} 241. 241. 1 4_light
        date
        echo ${a}: all files processed
    done
done


python ms_input.py poletti2010/exp2/cond7/12_8_pi_arc/1 exp2_cond7_12_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/14_8_pi_arc/1 exp2_cond7_14_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/15_8_pi_arc/1 exp2_cond7_15_8_pi_arc_1 1000 
wait
date
echo fifth three files processed