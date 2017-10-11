#!/bin/bash
#read video name

date
#python ms_input.py oof/detail detail 400 
python ms_network2.py oof/detail detail 240 40

if [ ! -d ./video/img_input/oof/cf_1 ]; then
        mkdir -p ./video/img_input/oof/cf_1
fi
if [ ! -d ./video/img_input/oof/cf_1 ]; then
        mkdir -p ./video/img_input/oof/cf_1
fi

python displacement_oof.py
python ms_input.py oof/detail detail 400 
python ms_network2.py oof/detail detail 240 40

date
echo done

for i in {1..10} #rereplace by 10
do
    if [ ! -d ./video/img_input/oof/cf_1 ]; then
            mkdir -p ./video/img_input/oof/cf_1
    fi
    if [ ! -d ./video/img_input/oof/cf_1 ]; then
            mkdir -p ./video/img_input/oof/cf_1
    fi 
    python ms_input.py oof/cf_${i} cf_${i} 700 
done

    python ms_input.py oof/cf_${i} cf_${i} 400 
    python ms_network2.py oof/cf_${i} cf_${i} 240 40