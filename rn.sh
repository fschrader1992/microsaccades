for k in 1 #2 3 4 5 6 7
do

    for v in 6 #0 3 6
    do
        for i in 30 #45 60 75 105 #0
        do
            for q in 64 #32 16 8 #4 2 96
            do
                #./data/jitter/chess_size${q}_vel${v}on105off${i}/${k}
                #rename 's/\off0/\off${i}/' data/jitter/chess_size8_vel6on105off${i}/${k}/*chess_size*_vel6on105off0*.data
                #for j in data/jitter/chess_size${q}_vel6on105off${i}/${k}/*chess_size*_vel6on105off0*.data; do
                #    new=$(printf "off%d" "$i") #04 pad to length of 4
                #    mv -i -- "$j" "$new"
                done
            done
        done
    done
done