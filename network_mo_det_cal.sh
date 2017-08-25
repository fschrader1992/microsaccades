#!/bin/bash
#read video name

for i in 1 10 #100
do
    echo '--- DELAY ${i} STARTS NOW ---'
    date
    python ms_network.py mo_det_cal mo_det_cal_1fr 121. 21. $i 1.
    #python ms_network.py mo_det_cal mo_det_cal_2fr 121. 21. $i 2.
    #python ms_network.py mo_det_cal mo_det_cal_3fr 121. 21. $i 3.
    #python ms_network.py mo_det_cal mo_det_cal_5fr 121. 21. $i 5.

    echo 'four files processed'

    #python ms_network.py mo_det_cal mo_det_cal_10fr 121. 21. $i 10.
    #python ms_network.py mo_det_cal mo_det_cal_16fr 121. 21. $i 16.
    #python ms_network.py mo_det_cal mo_det_cal_20fr 121. 21. $i 20.
    #python ms_network.py mo_det_cal mo_det_cal_33fr 121. 21. $i 33.

    echo '---FINISHED---'
done

date

python graphics_mo_det_cal.py
echo done