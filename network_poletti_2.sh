#!/bin/bash

date
echo 'STARTING SIMULATIONS'

for i in 8 9 10
do
    if [ ! -d ./data/poletti2010/exp2/cond3/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond3/network/${i};
    fi
    if [ ! -d ./data/poletti2010/exp1/cond1_light/network/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1_light/network/${i};
    fi

    date
    python ms_network_rest.py poletti2010/exp2/cond3/ ${i} ${i} exp2_cond3_${i} 241. 241. 2 3
    python ms_network.py poletti2010/exp1/cond1_light/ ${i} ${i} exp1_cond1_light_${i} 241. 241. 1 1_light
    date
    echo ${i}: seven files processed
done

for i in 2 3 4 5 6 7
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
    python ms_network.py poletti2010/exp1/cond3/ 1 ${i} exp1_cond3 241. 241. 1 3
    python ms_network.py poletti2010/exp2/cond5/ 1 ${i} exp2_cond5 241. 241. 2 5
    date
    echo ${i}: two files processed
done

echo 'exp 1 con 1, 1_light, 3 and 3 simulated'