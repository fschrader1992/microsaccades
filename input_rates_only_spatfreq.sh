#!/bin/bash
#read video name

date
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg1spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg2spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg5spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg10spat 1000 

wait
echo four files processed
date

python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg15spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg20spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg30spat 1000 &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg60spat 1000

wait
echo four files processed

#python graphics_spat_freq.py

#echo graphics created
date
echo done spat_freq only