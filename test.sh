#!/bin/sh

python ms_input.py spat_freq spatfreq_0fr0deg10spat & 
python ms_input.py spat_freq spatfreq_0fr0deg2spat

wait
echo all processes complete