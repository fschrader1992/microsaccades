#!/bin/bash
#read video name

date
for files in video/spat_freq/*; do
    filename=$(basename "$files")
    filename="${filename%.*}"
    echo $filename
    python ms_input.py spat_freq $filename
done
date
for files in video/phases/*; do
    filename=$(basename "$files")
    filename="${filename%.*}"
    echo $filename
    python ms_input.py phases $filename
done
date
echo "done"