#!/bin/bash
#read video name

for i in {0..330..30}
do
    deg=$i
    echo $deg
    if [ ! -d ./video/img_input/phases_0fr${deg}deg ]; then
    mkdir -p ./video/img_input/phases_0fr${deg}deg;
    fi
    if [ ! -d ./img/video/phases_0fr${deg}deg ]; then
    mkdir -p ./img/video/phases_0fr${deg}deg;
    fi
    echo "folders created"
    python image_creation_phase.py $deg
    echo "images created"
    cd ./video/img_input/phases_0fr${deg}deg
    ffmpeg -f image2 -r 24 -i second%3d.png -vcodec mpeg4 -y ../../phases/phases_0fr${deg}deg.mp4
    echo "simple video created"
    cd ../../..
done
echo "done"