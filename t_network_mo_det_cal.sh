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

#if [ ! -d ./video/img_input/mo_det_cal/mo_det_cal_dot_1fr ]; then
#    mkdir -p ./video/img_input/mo_det_cal/mo_det_cal_dot_1fr;
#fi
#if [ ! -d ./data/mo_det_cal/mo_det_cal_dot_1fr ]; then
#    mkdir -p ./data/mo_det_cal/mo_det_cal_dot_1fr;
#fi
#python image_creation_dot.py mo_det_cal/mo_det_cal_dot_1fr 1 0 0 0 0 
#python ms_input.py mo_det_cal/mo_det_cal_dot_0fr mo_det_cal_dot_0fr 400 

for i in 34 #10 100
do
    echo --- DELAY ${i} STARTS NOW ---
    date
    echo delay: $i
    python ms_network.py mo_det_cal/mo_det_cal_dot_0fr mo_det_cal_dot_0fr 121. 31. 34 0. 40.
    #echo delay: $i
    python ms_network.py mo_det_cal/mo_det_cal_dot_1fr mo_det_cal_dot_1fr 121. 31. 34 1. 40.
    echo delay: $i
    #python ms_network.py mo_det_cal/mo_det_cal_dot_m1fr mo_det_cal_dot_m1fr 121. 31. 34 -1. 40.
    #echo delay: $i
    python ms_network.py mo_det_cal/mo_det_cal_dot_16fr mo_det_cal_dot_16fr 121. 31. 34 16. 40.
    #echo delay: $i
    #python ms_network.py mo_det_cal/mo_det_cal_dot_m16fr mo_det_cal_dot_m16fr 121. 31. 134 -16. 40. 
    echo delay: $i
    python ms_network.py mo_det_cal/mo_det_cal_dot_33fr mo_det_cal_dot_33fr 121. 31. 34 33. 40.
    #echo delay: $i
    #python ms_network.py mo_det_cal/mo_det_cal_dot_m33fr mo_det_cal_dot_m33fr 121. 31. 34 -33. 40.
    
    #echo delay: $i
    #python ms_network.py mo_det_cal/mo_det_cal_10fr mo_det_cal_10fr 121. 21. 2 1. 
    #echo delay: $i
    #python ms_network.py mo_det_cal/mo_det_cal_66fr mo_det_cal_66fr 121. 21. 2 1. 

    echo ---
   
done

echo '---FINISHED---'
date

#python graphics_mo_det_cal.py
#echo graphics done