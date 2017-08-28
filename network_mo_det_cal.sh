#!/bin/bash
#read video name

for v in 1 2 3 5 10 16 20 33 
do

    if [ ! -d ./data/mo_det_cal/mo_det_cal_${v}fr ]; then
        mkdir -p ./data/mo_det_cal/mo_det_cal_${v}fr;
    fi
    python ms_input.py mo_det_cal/mo_det_cal_${v}fr mo_det_cal_${v}fr 300 

done

for v in 5 10 20 30 
do
    for i in 1 10 100
    do
        echo --- DELAY ${i} WITH WEIGHT ${v} STARTS NOW ---
        date
        python ms_network.py mo_det_cal/mo_det_cal_1fr mo_det_cal_1fr 121. 21. $i 1. $v
        python ms_network.py mo_det_cal_2fr mo_det_cal_2fr 121. 21. $i 2. $v
        python ms_network.py mo_det_cal_3fr mo_det_cal_3fr 121. 21. $i 3. $v
        python ms_network.py mo_det_cal_5fr mo_det_cal_5fr 121. 21. $i 5. $v

        echo delay ${i} weight ${v}: first four files processed

        python ms_network.py mo_det_cal_10fr mo_det_cal_10fr 121. 21. $i 10. $v
        python ms_network.py mo_det_cal_16fr mo_det_cal_16fr 121. 21. $i 16. $v
        python ms_network.py mo_det_cal_20fr mo_det_cal_20fr 121. 21. $i 20. $v
        python ms_network.py mo_det_cal_33fr mo_det_cal_33fr 121. 21. $i 33. $v
    
        echo delay ${i} weight ${v}: second four files processed
        echo ---
        
    done
done
echo '---FINISHED---'
date

python graphics_mo_det_cal.py
echo done