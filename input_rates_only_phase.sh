#!/bin/bash
#read video name

date
python ms_input.py phases phases_0fr0deg &
python ms_input.py phases phases_0fr30deg &
python ms_input.py phases phases_0fr60deg &
python ms_input.py phases phases_0fr90deg &
python ms_input.py phases phases_0fr120deg &
python ms_input.py phases phases_0fr150deg &

wait
echo six files processed

date
python ms_input.py phases phases_0fr180deg &
python ms_input.py phases phases_0fr210deg &
python ms_input.py phases phases_0fr240deg &
python ms_input.py phases phases_0fr270deg &
python ms_input.py phases phases_0fr300deg &
python ms_input.py phases phases_0fr330deg &

wait
echo six files processed
date
echo done