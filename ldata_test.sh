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

for a in 3 
do
    for i in 1 #{1..10} 
    do
        date
        
        #python together_new.py poletti2010/exp1/cond2/ ${a}_8_pi_arc ${i} exp1_cond2
        
        #python ldata_test.py poletti2010/exp2/cond1/ ${a}_8_pi_arc ${i} exp2_cond1 2 1
        
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


for i in 4
do
    date
    python ldata_test.py poletti2010/exp2/cond1/ ${i} exp2_cond1_${i} 2 1
    echo ${i}: seven files processed
done


echo '---FINISHED---'
date

#python graphics_mo_det_cal.py
#echo graphics done