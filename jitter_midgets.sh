for k in 1 #2 3 4 5 6 7
do
    for v in 6 #0 3 6
    do
        for i in 0 15 30 45 60 75 90 105 #0 15 30 45 60 75 90 105 
        do
            for q in 16 #64 32 16 8 #4 2 96 
            do
                for m in 4 #2 4 6 8 
                do
                echo chess_large: $v $q
                #python graphics_chess.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $q
                #python ms_network_no_m.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $m 200 #315
                #python ms_network.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $m 200 #315
                
                python ms_network.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i
                python together_lr.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $m
                
                if [ ! -d ./img/animation/jitter/chess_size${q}_vel${v}on105off${i} ]; then
                    mkdir -p ./img/animation/jitter/chess_size${q}_vel${v}on105off${i}
                fi
                python ms_network_midget_inspect.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 130. 130. $v $i $m #481. 121.
                python video_evaluation.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $m jitter
                
                done
            done
        done
    done
done