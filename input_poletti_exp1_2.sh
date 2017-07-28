#!/bin/bash

#experiment1 cond3 and exp2 cond5 only need one video
python image_creation_dot.py exp1/cond3 0 0 0 0 0 
python image_creation_dot.py exp2/cond5 0 0 0 1 0 
echo '2 and 5 created'

#exp1 cond1 dark, exp1 cond 1 light and exp1 cond 3 light has just eye movements
for i in {1..10} #rereplace by 10
do
    #exp1
    if [ ! -d ./video/img_input/poletti2010/exp1/cond1/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond1/${i};
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond1_light/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond1_light/${i};
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond3_light/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond3_light/${i};
    fi
    python image_creation_dot.py exp1/cond1/${i} 0 0 1 0 0 
    python image_creation_dot_light.py exp1/cond1_light/${i} 0 0 1 0 1 
    python image_creation_dot_light.py exp1/cond3_light/${i} 0 0 0 0 1 
    
    if [ ! -d ./data/poletti2010/exp1/cond1/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1/${i};
    fi
    if [ ! -d ./data/poletti2010/exp1/cond1_light/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1_light/${i};
    fi
    if [ ! -d ./data/poletti2010/exp1/cond3_light/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond3_light/${i};
    fi
    
    #exp2
    if [ ! -d ./video/img_input/poletti2010/exp2/cond1/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond1/${i};
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond3/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond3/${i};
    fi
    python image_creation_dot.py exp2/cond1/${i} 0 0 1 1 1 
    python image_creation_dot.py exp2/cond3/${i} 0 0 0 1 1 
    
    if [ ! -d ./data/poletti2010/exp2/cond1/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond1/${i};
    fi
    if [ ! -d ./data/poletti2010/exp2/cond3/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond3/${i};
    fi
done
 
 
for a in {0..15} #15
do
    date
    if [ ! -d ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc;
    fi
    
    if [ ! -d ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc;
    fi
    
    
    
    if [ ! -d ./data/poletti2010/exp1/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc;
    fi
    
    if [ ! -d ./data/poletti2010/exp2/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp2/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond4/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp2/cond6/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond6/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp2/cond7/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond7/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp2/cond8/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond8/${a}_8_pi_arc;
    fi
    echo "folders created"
    
    #exp1 cond4 is 
    python image_creation_dot.py exp1/cond4/${a}_8_pi_arc 1 $a 0 0 0
    
    for i in 1 #{1..10} 
    do
        #exp1
        if [ ! -d ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc/${i};
        fi
        python image_creation_dot.py exp1/cond2/${a}_8_pi_arc/${i} 1 $a 1 0 0 
        python image_creation_dot_light.py exp1/cond2_light/${a}_8_pi_arc/${i} 1 $a 1 0 1
        python image_creation_dot_light.py exp1/cond4_light/${a}_8_pi_arc/${i} 1 $a 0 0 1
        
        #exp2
        if [ ! -d ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc/${i};
        fi
        python image_creation_dot.py exp2/cond2/${a}_8_pi_arc/${i} 1 $a 1 1 1 
        python image_creation_dot.py exp2/cond4/${a}_8_pi_arc/${i} 1 $a 0 1 1
        python image_creation_dot.py exp2/cond6/${a}_8_pi_arc/${i} 1 $a 0 1 0
        python image_creation_dot_7_8.py exp2/cond6/${a}_8_pi_arc/${i} 0 $a 1 1 1 1
        python image_creation_dot_7_8.py exp2/cond6/${a}_8_pi_arc/${i} 1 $a 1 1 1 1
        
    done
    
    echo ${a}: images created

done

date
echo 'STARTING SIMULATIONS'

python ms_input.py poletti2010/exp1/cond3 exp1_cond3 1000 &
python ms_input.py poletti2010/exp2/cond5 exp2_cond5 1000 &
python ms_input.py poletti2010/exp1/cond4/0_8_pi_arc exp1_cond4_0_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/1_8_pi_arc exp1_cond4_1_8_pi_arc 1000 

wait

echo '2, 4 and 5 part 1 simulated'

python ms_input.py poletti2010/exp1/cond4/2_8_pi_arc exp1_cond4_2_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/3_8_pi_arc exp1_cond4_3_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/14_8_pi_arc exp1_cond4_14_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/15_8_pi_arc exp1_cond4_15_8_pi_arc 1000 

wait

echo '2, 4 and 5 completely simulated'
date

python ms_input.py poletti2010/exp1/cond4/5_8_pi_arc exp1_cond4_5_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/7_8_pi_arc exp1_cond4_7_8_pi_arc 1000 &
python ms_input.py poletti2010/exp1/cond4/10_8_pi_arc exp1_cond4_10_8_pi_arc 1000 & 
python ms_input.py poletti2010/exp1/cond4/12_8_pi_arc exp1_cond4_12_8_pi_arc 1000 &

wait

echo 'condition 4 also simulated'

for a in 0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_input.py poletti2010/exp1/cond2/${a}_8_pi_arc/${i} exp1_cond2_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp1/cond2_light/${a}_8_pi_arc/${i} exp1_cond2_light_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp1/cond4_light/${a}_8_pi_arc/${i} exp1_cond4_light_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond2/${a}_8_pi_arc/${i} exp2_cond2_${a}_8_pi_arc_${i} 1000 
        wait
        date
        echo ${a}: first four files processed
        
        python ms_input.py poletti2010/exp2/cond4/${a}_8_pi_arc/${i} exp2_cond4_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond6/${a}_8_pi_arc/${i} exp2_cond6_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond7/${a}_8_pi_arc/${i} exp2_cond7_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond8/${a}_8_pi_arc/${i} exp2_cond8_${a}_8_pi_arc_${i} 1000 
        wait
        date
        echo ${a}: second four files processed
    done
done

for i in {1..10} 
do
    date
    python ms_input.py poletti2010/exp1/cond1/${i} exp1_cond1_${i} 1000 &
    python ms_input.py poletti2010/exp1/cond1_light/${i} exp1_cond1_light_${i} 1000 &
    python ms_input.py poletti2010/exp1/cond3_light/${i} exp1_cond3_light_${i} 1000 &
    python ms_input.py poletti2010/exp2/cond1/${i} exp2_cond1_${i} 1000 &
    python ms_input.py poletti2010/exp2/cond3/${i} exp2_cond3_${i} 1000 
    wait
    date
    echo ${i}: five files processed
done

echo 'exp 1 con 1, 1_light and 3 simulated'

date
echo done with this