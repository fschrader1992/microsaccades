#!/bin/bash
#read video name

date

#python image_creation_spatfreq.py 5
#echo 10
#python image_creation_spatfreq.py 10
#echo 15
#python image_creation_spatfreq.py 15
#echo 20
#python image_creation_spatfreq.py 20
#echo 30
#python image_creation_spatfreq.py 30
#echo 60
#python image_creation_spatfreq.py 60

#echo image input created

python ms_input.py phases/phases_0fr0deg_32px_only_border phases_0fr0deg_32px_only_border 1000

echo phases done

date
python ms_input.py spat_freq/spatfreq_0fr0deg1spat 1spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg2spat 2spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg5spat 5spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg10spat 10spat 400

wait
echo four files processed

date
python ms_input.py spat_freq/spatfreq_0fr0deg15spat 15spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg20spat 20spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg30spat 30spat 400 &
python ms_input.py spat_freq/spatfreq_0fr0deg60spat 60spat 400

wait
echo four files processed

#python graphics_spat_freq.py

#echo graphics created
date
echo done