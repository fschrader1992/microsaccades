#!/bin/bash
#read video name

#for v in 16 #  16 33 
#do

#   if [ ! -d ./data/mo_det_cal/mo_det_cal_dot_${v}fr ]; then
#       mkdir -p ./data/mo_det_cal/mo_det_cal_dot_${v}fr;
#   fi
#   python ms_input.py mo_det_cal/mo_det_cal_dot_${v}fr mo_det_cal_dot_${v}fr 400 
   
#   if [ ! -d ./data/mo_det_cal/mo_det_cal_dot_m${v}fr ]; then
#       mkdir -p ./data/mo_det_cal/mo_det_cal_dot_m${v}fr;
#   fi
#   python ms_input.py mo_det_cal/mo_det_cal_dot_m${v}fr mo_det_cal_dot_m${v}fr 400 

#done

#if [ ! -d ./video/img_input/mo_det_cal/mo_det_cal_dot_deg60_16fr ]; then
#    mkdir -p ./video/img_input/mo_det_cal/mo_det_cal_dot_deg60_16fr;
#fi
#if [ ! -d ./data/mo_det_cal/mo_det_cal_dot_deg60_16fr ]; then
#    mkdir -p ./data/mo_det_cal/mo_det_cal_dot_deg60_16fr;
#fi
#python image_creation_dot.py mo_det_cal/mo_det_cal_dot_deg60_16fr 1 0 0 0 0 
#python ms_input.py mo_det_cal/mo_det_cal_dot_deg60_16fr mo_det_cal_dot_deg60_16fr 400 

#for i in 1 #50. 100. 200.  #10 100
#do

for a in 1 5 7 10 12 15  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        
        #python together_cells.py poletti2010/exp1/cond2/ ${a}_8_pi_arc ${i} exp1_cond2
        
        python together_cells.py poletti2010/exp1/cond2/ ${a}_8_pi_arc ${i} exp1_cond2 1 2
        python together_cells.py poletti2010/exp1/cond2_light/ ${a}_8_pi_arc ${i} exp1_cond2_light 1 2_light
        python together_cells.py poletti2010/exp1/cond4/ ${a}_8_pi_arc ${i} exp1_cond4 1 4
        python together_cells.py poletti2010/exp2/cond2/ ${a}_8_pi_arc ${i} exp2_cond2 2 2
        python together_cells.py poletti2010/exp2/cond4/ ${a}_8_pi_arc ${i} exp2_cond4 2 4
        python together_cells.py poletti2010/exp2/cond6/ ${a}_8_pi_arc ${i} exp2_cond6 2 6
        python together_cells.py poletti2010/exp2/cond8/ ${a}_8_pi_arc ${i} exp2_cond8 2 8
        
        echo --- Simulation ${i} STARTS NOW ---
        date
        echo number: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_rem_0fr mo_det_cal_dot_rem_0fr 61. 16. 34 0. 40. 
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_1fr mo_det_cal_dot_1fr 121. 31. 34 1. 40. 
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_m1fr mo_det_cal_dot_m1fr 121. 31. 34 -1. 40.
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_16fr 1 ${i} mo_det_cal_dot_16fr 121. 31. 1 1
        #python ms_network.py mo_det_cal/mo_det_cal_dot_deg90_16fr mo_det_cal_dot_deg90_16fr 31. 121. 34 16. 40. 378.
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_m16fr mo_det_cal_dot_m16fr 121. 31. 134 -16. 40. 
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_33fr mo_det_cal_dot_33fr 121. 31. 34 33. 40. 379.
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_m33fr mo_det_cal_dot_m33fr 121. 31. 34 -33. 40.
        
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_1fr mo_det_cal_dot_1fr 121. 31. 34 16. 40. 
        #python ms_network.py mo_det_cal/mo_det_cal_dot_m16fr mo_det_cal_dot_m16fr 121. 31. 34 16. 40. 
        #python ms_network.py mo_det_cal/mo_det_cal_dot_16fr mo_det_cal_dot_16fr 121. 31. 34 16. 40. 
        #echo delay: $i
        #python ms_network.py mo_det_cal/mo_det_cal_dot_33fr mo_det_cal_dot_33fr 121. 31. 34 33. 40.

    echo ---
    done
   
done

for a in 1 2 3 5 7 10  #0 1 2 3 5 7 10 12 14 15 
do
    for i in 1 #{1..10} 
    do
        date
        
        python together_cells.py poletti2010/exp1/cond4_light/ ${a}_8_pi_arc ${i} exp1_cond4_light 1 4_light
        python together_cells.py poletti2010/exp2/cond7/ ${a}_8_pi_arc ${i} exp2_cond7 2 7

    echo ---
    done
   
done

for i in 1 2 3 4 5 6 7 8 9 10 
do
    date
    python together_cells.py poletti2010/exp1/cond1/ ${i} exp1_cond1_${i} 1 1
    python together_cells.py poletti2010/exp1/cond1_light/ ${i} exp1_cond1_light_${i} 1 1_light
    python together_cells.py poletti2010/exp1/cond3/ ${i} exp1_cond3 1 3
    python together_cells.py poletti2010/exp1/cond3_light/ ${i} exp1_cond3_light_${i} 1 3_light
    python together_cells.py poletti2010/exp2/cond1/ ${i} exp2_cond1_${i} 2 1
    python together_cells.py poletti2010/exp2/cond3/ ${i} exp2_cond3_${i} 2 3
    python together_cells.py poletti2010/exp2/cond5/ ${i} exp2_cond5 2 5
    date
    echo ${i}: seven files processed
done


echo '---FINISHED---'
date

#python graphics_mo_det_cal.py
#echo graphics done