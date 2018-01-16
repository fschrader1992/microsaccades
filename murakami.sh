#!/bin/bash

for k in 8 9 10 11 #2 3 4 5 6 7
do

    for i in 0 15 30 45 60 75 90 105
    do
        echo off duty ${i}ms
        if [ ! -d ./video/img_input/murakami/on105off${i}/${k} ]; then
            mkdir -p ./video/img_input/murakami/on105off${i}/${k}
        fi
        if [ ! -d ./data/murakami/on105off${i}/${k} ]; then
            mkdir -p ./data/murakami/on105off${i}/${k}
        fi
        if [ ! -d ./data/murakami/on105off${i}/network/${k} ]; then
            mkdir -p ./data/murakami/on105off${i}/network/${k};
        fi
        python image_creation_murakami.py on105off${i}/${k} 105 $i 500
    done

    echo creation done 

    python ms_input.py murakami/on105off0/${k} on105off0 500 &
    python ms_input.py murakami/on105off15/${k} on105off15 500 &
    python ms_input.py murakami/on105off30/${k} on105off30 500 &
    python ms_input.py murakami/on105off45/${k} on105off45 500 &
    wait

    python ms_input.py murakami/on105off60/${k} on105off60 500 &
    python ms_input.py murakami/on105off75/${k} on105off75 500 &
    python ms_input.py murakami/on105off90/${k} on105off90 500 &
    python ms_input.py murakami/on105off105/${k} on105off105 500 &
    wait

    echo input done
        
    for i in 0 15 30 45 60 75 90 105
    do
        python ms_network_no_m.py murakami/on105off${i} ${k} on105off${i} 481. 121. 105 $i
        python ms_network_rest.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        python ms_network_rest_modet.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
    done
done
