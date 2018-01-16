for k in 2 3 4 5 6 7
do
    for i in 1 2 3 4 5 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40
    do
        echo Grating ${i}
        for j in 3 6
        do
            echo Grating ${i}, velocity ${j}
            if [ ! -d ./data/mtesting/mt_v${j}_c${i}/network/30/${k} ]; then
                mkdir -p ./data/mtesting/mt_v${j}_c${i}/network/30/${k};
            fi
            if [ ! -d ./data/mtesting/mt_v${j}_c${i}/network/15/${k} ]; then
                mkdir -p ./data/mtesting/mt_v${j}_c${i}/network/15/${k};
            fi
        done
        python ms_network.py mtesting/mt_v3_c${i} 30/${k} mt_v3_c${i} 121. 31. $i 3 30
        python ms_network.py mtesting/mt_v6_c${i} 30/${k} mt_v6_c${i} 121. 31. $i 6 30
        python ms_network.py mtesting/mt_v6_c${i} 15/${k} mt_v6_c${i} 121. 31. $i 6 15
    done
done
