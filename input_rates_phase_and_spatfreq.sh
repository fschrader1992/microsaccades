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
echo eightfiles processed

date
python ms_input.py phases phases_0fr0deg_24px &
python ms_input.py phases phases_0fr0deg_24px_only_border &
python ms_input.py phases phases_0fr0deg &
python ms_input.py phases phases_0fr30deg &
python ms_input.py phases phases_0fr60deg &
python ms_input.py phases phases_0fr90deg 

wait
echo six files processed

date
python ms_input.py phases phases_0fr120deg &
python ms_input.py phases phases_0fr150deg &
python ms_input.py phases phases_0fr180deg &
python ms_input.py phases phases_0fr210deg &
python ms_input.py phases phases_0fr240deg &
python ms_input.py phases phases_0fr270deg &
python ms_input.py phases phases_0fr300deg &
python ms_input.py phases phases_0fr330deg 

wait
echo six files processed
date
echo done