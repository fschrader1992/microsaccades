for i in 1 2 3 4 5 8 10 15 20
do
    echo Grating ${i}
    for j in 6 
    do
        echo Grating ${i}, velocity ${j}
        if [ ! -d ./video/img_input/mtesting/mt_v${j}_c${i} ]; then
            mkdir -p ./video/img_input/mtesting/mt_v${j}_c${i}
        fi
        if [ ! -d ./data/mtesting/mt_v${j}_c${i} ]; then
            mkdir -p ./data/mtesting/mt_v${j}_c${i}
        fi
        if [ ! -d ./data/mtesting/mt_v${j}_c${i}/network/30 ]; then
            mkdir -p ./data/mtesting/mt_v${j}_c${i}/network/30;
        fi
        python image_creation_midget_testing.py mt_v${j}_c${i} $i $j
    done
    python ms_input.py mtesting/mt_v3_c${i} mt_v3_c${i} 400 &
    python ms_input.py mtesting/mt_v6_c${i} mt_v6_c${i} 400 
    wait
    python ms_network.py mtesting/mt_v3_c${i} 30 mt_v3_c${i} 121. 31. $i 3 30
    python ms_network.py mtesting/mt_v6_c${i} 30 mt_v6_c${i} 121. 31. $i 6 30
done

for i in 1 2 3 4 5 8 10 15 20
do
    echo Grating ${i}
    for j in 12 18 
    do
        if [ ! -d ./video/img_input/mtesting/mt_v${j}_c${i} ]; then
            mkdir -p ./video/img_input/mtesting/mt_v${j}_c${i}
        fi
        if [ ! -d ./data/mtesting/mt_v${j}_c${i} ]; then
            mkdir -p ./data/mtesting/mt_v${j}_c${i}
        fi
        if [ ! -d ./data/mtesting/mt_v${j}_c${i}/network/15 ]; then
            mkdir -p ./data/mtesting/mt_v${j}_c${i}/network/15;
        fi
        echo Grating ${i}, velocity 0.01*${j}
        python image_creation_midget_testing.py mt_v${j}_c${i} $i $j
    done
    python ms_input.py mtesting/mt_v12_c${i} mt_v12_c${i} 400 &
    python ms_input.py mtesting/mt_v18_c${i} mt_v18_c${i} 400 
    wait
    if [ ! -d ./data/mtesting/mt_v6_c${i}/network/15 ]; then
        mkdir -p ./data/mtesting/mt_v6_c${i}/network/15;
    fi
    python ms_network.py mtesting/mt_v6_c${i} 15 mt_v6_c${i} 121. 31. $i 6 15
    python ms_network.py mtesting/mt_v12_c${i} 15 mt_v12_c${i} 121. 31. $i 12 15
    python ms_network.py mtesting/mt_v18_c${i} 15 mt_v18_c${i} 121. 31. $i 18 15
done