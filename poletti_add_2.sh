#!/bin/bash

#experiment1 cond3 and exp2 cond5 only need one video
python image_creation_dot.py exp1/cond3 0 0 0 0 0 
python image_creation_dot.py exp2/cond5 0 0 0 1 0 
echo '2 and 5 created'

#exp1 cond1 dark, exp1 cond 1 light and exp1 cond 3 light has just eye movements
for i in 8 9 10 11 #rereplace by 10
do
    #exp1
    if [ ! -d ./video/img_input/poletti2010/exp1/cond1/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond1/${i};
    fi
    #if [ ! -d ./video/img_input/poletti2010/exp1/cond1_light/${i} ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp1/cond1_light/${i};
    #fi
    #if [ ! -d ./video/img_input/poletti2010/exp1/cond3_light/${i} ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp1/cond3_light/${i};
    #fi
    python image_creation_dot.py exp1/cond1/${i} 0 0 1 0 0 
    #python image_creation_dot_light.py exp1/cond1_light/${i} 0 0 1 0 1 
    #python image_creation_dot_light.py exp1/cond3_light/${i} 0 0 0 0 1 
    
    if [ ! -d ./data/poletti2010/exp1/cond1/${i} ]; then
        mkdir -p ./data/poletti2010/exp1/cond1/${i};
    fi
    #if [ ! -d ./data/poletti2010/exp1/cond1_light/${i} ]; then
    #    mkdir -p ./data/poletti2010/exp1/cond1_light/${i};
    #fi
    #if [ ! -d ./data/poletti2010/exp1/cond3_light/${i} ]; then
    #    mkdir -p ./data/poletti2010/exp1/cond3_light/${i};
    #fi
    
    #exp2
    if [ ! -d ./video/img_input/poletti2010/exp2/cond1/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond1/${i};
    fi
    #if [ ! -d ./video/img_input/poletti2010/exp2/cond3/${i} ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp2/cond3/${i};
    #fi
    python image_creation_dot.py exp2/cond1/${i} 0 0 1 1 1 
    #python image_creation_dot.py exp2/cond3/${i} 0 0 0 1 1 
    
    if [ ! -d ./data/poletti2010/exp2/cond1/${i} ]; then
        mkdir -p ./data/poletti2010/exp2/cond1/${i};
    fi
    #if [ ! -d ./data/poletti2010/exp2/cond3/${i} ]; then
    #    mkdir -p ./data/poletti2010/exp2/cond3/${i};
    #fi
done
 
 
for a in 0 2 8 13 #15
do
    date
    if [ ! -d ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc;
    fi
    #if [ ! -d ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc;
    #fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    #if [ ! -d ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc;
    #fi
    
    if [ ! -d ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc;
    fi
    #if [ ! -d ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc;
    #fi
    #if [ ! -d ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc;
    #fi
    #if [ ! -d ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc ]; then
    #    mkdir -p ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc;
    #fi
    
    
    
    if [ ! -d ./data/poletti2010/exp1/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond2/${a}_8_pi_arc;
    fi
    #if [ ! -d ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc ]; then
    #    mkdir -p ./data/poletti2010/exp1/cond2_light/${a}_8_pi_arc;
    #fi
    if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    #if [ ! -d ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc ]; then
    #    mkdir -p ./data/poletti2010/exp1/cond4_light/${a}_8_pi_arc;
    #fi
    
    if [ ! -d ./data/poletti2010/exp2/cond2/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond2/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp2/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp2/cond4/${a}_8_pi_arc;
    fi
    #if [ ! -d ./data/poletti2010/exp2/cond6/${a}_8_pi_arc ]; then
    #    mkdir -p ./data/poletti2010/exp2/cond6/${a}_8_pi_arc;
    #fi
    #if [ ! -d ./data/poletti2010/exp2/cond7/${a}_8_pi_arc ]; then
    #    mkdir -p ./data/poletti2010/exp2/cond7/${a}_8_pi_arc;
    #fi
    #if [ ! -d ./data/poletti2010/exp2/cond8/${a}_8_pi_arc ]; then
    #    mkdir -p ./data/poletti2010/exp2/cond8/${a}_8_pi_arc;
    #fi
    echo "folders created"
    
    #exp1 cond4 is 
    #python image_creation_dot.py exp1/cond4/${a}_8_pi_arc 1 $a 0 0 0
    
    for i in 1 #{1..10} 
    do
        #exp1
        if [ ! -d ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond2/${a}_8_pi_arc/${i};
        fi
        #if [ ! -d ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc/${i} ]; then
        #    mkdir -p ./video/img_input/poletti2010/exp1/cond2_light/${a}_8_pi_arc/${i};
        #fi
        if [ ! -d ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc/${i};
        fi
        #if [ ! -d ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc/${i} ]; then
        #    mkdir -p ./video/img_input/poletti2010/exp1/cond4_light/${a}_8_pi_arc/${i};
        #fi
        python image_creation_dot.py exp1/cond2/${a}_8_pi_arc/${i} 1 $a 1 0 0 
        #python image_creation_dot_light.py exp1/cond2_light/${a}_8_pi_arc/${i} 1 $a 1 0 1
        python image_creation_dot.py exp1/cond4/${a}_8_pi_arc/${i} 1 $a 0 0 0
        #python image_creation_dot_light.py exp1/cond4_light/${a}_8_pi_arc/${i} 1 $a 0 0 1
        
        #exp2
        if [ ! -d ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond2/${a}_8_pi_arc/${i};
        fi
        if [ ! -d ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp2/cond4/${a}_8_pi_arc/${i};
        fi
        #if [ ! -d ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc/${i} ]; then
        #    mkdir -p ./video/img_input/poletti2010/exp2/cond6/${a}_8_pi_arc/${i};
        #fi
        #if [ ! -d ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc/${i} ]; then
        #    mkdir -p ./video/img_input/poletti2010/exp2/cond7/${a}_8_pi_arc/${i};
        #fi
        #if [ ! -d ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc/${i} ]; then
        #    mkdir -p ./video/img_input/poletti2010/exp2/cond8/${a}_8_pi_arc/${i};
        #fi
        python image_creation_dot.py exp2/cond2/${a}_8_pi_arc/${i} 1 $a 1 1 1 
        python image_creation_dot.py exp2/cond4/${a}_8_pi_arc/${i} 1 $a 0 1 1
        #python image_creation_dot.py exp2/cond6/${a}_8_pi_arc/${i} 1 $a 0 1 0
        #python image_creation_dot_7_8.py exp2/cond7/${a}_8_pi_arc/${i} 0 $a 1 1 1 1
        #python image_creation_dot_7_8.py exp2/cond8/${a}_8_pi_arc/${i} 1 $a 1 1 1 1
        
    done
    
    echo ${a}: images created

done

for a in 1 5 7 12 #15
do
    date
    if [ ! -d ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    if [ ! -d ./data/poletti2010/exp1/cond4/${a}_8_pi_arc ]; then
        mkdir -p ./data/poletti2010/exp1/cond4/${a}_8_pi_arc;
    fi
    
    for i in 1 #{1..10} 
    do
        if [ ! -d ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc/${i} ]; then
            mkdir -p ./video/img_input/poletti2010/exp1/cond4/${a}_8_pi_arc/${i};
        fi
        python image_creation_dot.py exp1/cond4/${a}_8_pi_arc/${i} 1 $a 0 0 0
    done
    
    echo ${a}: images created

done

date
echo 'STARTING SIMULATIONS'

for a in 0 2 8 13 #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        python ms_input.py poletti2010/exp1/cond2/${a}_8_pi_arc/${i} exp1_cond2_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp1/cond4/${a}_8_pi_arc/${i} exp1_cond4_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond2/${a}_8_pi_arc/${i} exp2_cond2_${a}_8_pi_arc_${i} 1000 &
        python ms_input.py poletti2010/exp2/cond4/${a}_8_pi_arc/${i} exp2_cond4_${a}_8_pi_arc_${i} 1000 &
        wait
        date
        echo ${a}: four files processed
        
    done
done

python ms_input.py poletti2010/exp1/cond4/1_8_pi_arc/1 exp1_cond4_1_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4/5_8_pi_arc/1 exp1_cond4_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4/7_8_pi_arc/1 exp1_cond4_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4/12_8_pi_arc/1 exp1_cond4_12_8_pi_arc_1 1000 
wait

python ms_input.py poletti2010/exp1/cond1/11 exp1_cond1_11 1000 &
python ms_input.py poletti2010/exp2/cond1/11 exp2_cond1_11 1000 &
python ms_input.py poletti2010/exp1/cond1/8 exp1_cond1_8 1000 &
python ms_input.py poletti2010/exp2/cond1/8 exp2_cond1_8 1000 
wait

python ms_input.py poletti2010/exp1/cond1/9 exp1_cond1_9 1000 &
python ms_input.py poletti2010/exp2/cond1/9 exp2_cond1_9 1000 &
python ms_input.py poletti2010/exp1/cond1/10 exp1_cond1_10 1000 &
python ms_input.py poletti2010/exp2/cond1/10 exp2_cond1_10 1000 
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
