#!/bin/bash
#read video name

date
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg1spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg2spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg5spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg10spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg15spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg20spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg30spat &
python ms_input_spatfreq.py spat_freq spatfreq_0fr0deg60spat

wait
echo eight files processed

python graphics_spat_freq.py

echo graphics created
date
echo done