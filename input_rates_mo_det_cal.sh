#!/bin/bash
#read video name

date
python ms_input.py mo_det_cal/mo_det_cal_0fr mo_det_cal_0fr 300 &
python ms_input.py mo_det_cal/mo_det_cal_66fr mo_det_cal_66fr 300

#date
#python ms_input.py mo_det_cal/mo_det_cal_1fr mo_det_cal_1fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_2fr mo_det_cal_2fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_3fr mo_det_cal_3fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_5fr mo_det_cal_5fr 300 

wait
echo four files processed

#date
#python ms_input.py mo_det_cal/mo_det_cal_10fr mo_det_cal_10fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_16fr mo_det_cal_16fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_20fr mo_det_cal_20fr 300 &
#python ms_input.py mo_det_cal/mo_det_cal_33fr mo_det_cal_33fr 300 

#wait
#echo four files processed

date
echo done