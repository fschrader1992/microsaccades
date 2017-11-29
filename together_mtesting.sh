for k in 2 3 4 5 6 7
do
    for i in 1 2 3 4 5 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 40
    do
        python together_lr.py mt_v3_c${i} 30/${k} 3 ${i} 30
        python together_lr.py mt_v6_c${i} 30/${k} 6 ${i} 30
        python together_lr.py mt_v6_c${i} 15/${k} 6 ${i} 15
        #python together_lr.py mt_v12_c${i} 15 12 ${i}
        #python together_lr.py mt_v18_c${i} 15 18 ${i}
    done
done
