#!/bin/bash
#read video name

for i in .01 .1 1 10
do
    date
    python ms_network.py mo_det_cal mo_det_cal_1fr $i &
    python ms_network.py mo_det_cal mo_det_cal_2fr $i &
    python ms_network.py mo_det_cal mo_det_cal_3fr $i &
    python ms_network.py mo_det_cal mo_det_cal_5fr $i 

    wait
    echo four files processed

    date
    python ms_network.py mo_det_cal mo_det_cal_0_5fr $i &
    python ms_network.py mo_det_cal mo_det_cal_0_2fr $i &
    python ms_network.py mo_det_cal mo_det_cal_0_1fr $i &
    python ms_network.py mo_det_cal mo_det_cal_10fr $i &
    python ms_network.py mo_det_cal mo_det_cal_20fr $i 

    wait
    echo five files processed
done

date
echo done