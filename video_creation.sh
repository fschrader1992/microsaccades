#!/bin/bash
#read video name
echo "Hello, "$USER".  Video name:"
read video

if [ ! -d ./video/img_input/$video ]; then
  mkdir -p ./video/img_input/$video;
fi
if [ ! -d ./img/video/$video ]; then
  mkdir -p ./img/video/$video;
fi
echo "folders created"
python image_creation.py
echo "images created"
cd ./video/img_input/$video
ffmpeg -f image2 -r 24 -i second%3d.png -vcodec mpeg4 -y toloop.mp4
echo "simple video created"
for i in {1..4}; do printf "file '%s'\n" toloop.mp4 >> list.txt; done
ffmpeg -f concat -i list.txt -c copy ../../$video.mp4
cd ../../..
echo "done"