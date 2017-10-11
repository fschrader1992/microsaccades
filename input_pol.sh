#!/bin/bash

for i in {1..10} #rereplace by 10
do
    #exp1
    if [ ! -d ./video/img_input/poletti2010/exp1/cond1/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond1/${i};
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond1_light/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond1_light/${i};
    fi
    if [ ! -d ./video/img_input/poletti2010/exp1/cond3_light/${i} ]; then
        mkdir -p ./video/img_input/poletti2010/exp1/cond3_light/${i};
    fi
    python image_creation_dot.py exp1/cond1/${i} 0 0 1 0 0 
    python image_creation_dot_light.py exp1/cond1_light/${i} 0 0 1 0 1 
    python image_creation_dot_light.py exp1/cond3_light/${i} 0 0 0 0 1 