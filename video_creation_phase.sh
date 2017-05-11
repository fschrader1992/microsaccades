#!/bin/bash
#read video name

for i in {0..0}
do
    deg=$i
    echo $deg
    if [ ! -d ./video/img_input/phases_0fr${deg}deg_36px_only_border ]; then
    mkdir -p ./video/img_input/phases_0fr${deg}deg_36px_only_border;
    fi
    echo "folders created"
    python image_creation_phase.py $deg
    echo "images created"
    cd ./video/img_input/phases_0fr${deg}deg_36px_only_border
    ffmpeg -f image2 -r 24 -i second%3d.png -vcodec mpeg4 -y ../../phases/phases_0fr${deg}deg_36px_only_border.mp4
    echo "simple video created"
    cd ../../..
done
echo "done"