#!/bin/bash
#read video name

date
python ms_input.py spat_freq spatfreq_0fr0deg1spat &
python ms_input.py spat_freq spatfreq_0fr0deg2spat &
python ms_input.py spat_freq spatfreq_0fr0deg5spat &
python ms_input.py spat_freq spatfreq_0fr0deg10spat

wait
echo four files processed

date
python ms_input.py spat_freq spatfreq_0fr0deg15spat &
python ms_input.py spat_freq spatfreq_0fr0deg20spat &
python ms_input.py spat_freq spatfreq_0fr0deg30spat &
python ms_input.py spat_freq spatfreq_0fr0deg60spat

wait
echo four files processed

date
python ms_input.py phases phases_0fr0deg &
python ms_input.py phases phases_0fr30deg &
python ms_input.py phases phases_0fr60deg &
python ms_input.py phases phases_0fr90deg &

wait
echo four files processed
date
python ms_input.py phases phases_0fr120deg &
python ms_input.py phases phases_0fr150deg &
python ms_input.py phases phases_0fr180deg &
python ms_input.py phases phases_0fr210deg &

wait
echo four files processed
date
python ms_input.py phases phases_0fr240deg &
python ms_input.py phases phases_0fr270deg &
python ms_input.py phases phases_0fr300deg &
python ms_input.py phases phases_0fr330deg &

wait
echo four files processed
date
echo done