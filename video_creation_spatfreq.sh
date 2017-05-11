#!/bin/bash
#read video name


if [ ! -d ./img/video/spatfreq ]; then
mkdir -p ./img/video/spatfreq;
fi
for i in 1 2 5 10 15 20 30 60
do
    sf=$i
    echo $sf
    if [ ! -d ./video/img_input/spatfreq_0fr0deg${sf}spat ]; then
    mkdir -p ./video/img_input/spatfreq_0fr0deg${sf}spat;
    fi
    echo "folders created"
    #python image_creation_spatfreq.py $sf
    echo "images created"
    cd ./video/img_input/spatfreq_0fr0deg${sf}spat
    ffmpeg -f image2 -r 24 -i second%3d.png -vcodec mpeg4 -vf scale=240:240 -y ../../spat_freq/spatfreq_0fr0deg${sf}spat.mp4
    echo "simple video created"
    cd ../../..
done
echo "done"