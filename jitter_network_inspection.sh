for k in 5 #10 #2 3 4 5 6 7 8 9 10 11 #sim number
do
    for i in 15 #105 #30 60 75 90 105 #0 15 30 45 60 75 90 105 #off
    do
        for q in 16 #30 60 75 90 105 #0 15 30 45 60 75 90 105
        do
            for m in 4 #30 60 75 90 105 #0 15 30 45 60 75 90 105
            do
                python data_clean.py murakami/on105off${i}/ ${k} on105off${i}
                #python ms_network_rest.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
                #python ms_network_rest_modet.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
                #python ms_network_rest_weight.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
                #python ms_network_rest_modet_weight.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
                #python together_murakami.py murakami/on105off${i}/ ${k} on105off${i}
                #python together_murakami_modet.py murakami/on105off${i}/ ${k} on105off${i}
                #python together_murakami_weight.py murakami/on105off${i}/ ${k} on105off${i}
                #python together_murakami_modet_weight.py murakami/on105off${i}/ ${k} on105off${i}
                
                #python ms_network_rest_no_gmn.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
                #python together_murakami_no_gmn.py murakami/on105off${i}/ ${k} on105off${i}
                
        #        python ms_network_inspect.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i

                #python ms_network_inspect.py mo_det_cal/mo_det_cal_10fr ${k} mo_det_cal_10fr 131. 131. 105 $i
                #python video_evaluation.py mo_det_cal/mo_det_cal_10fr ${k} mo_det_cal_10fr 131. 131. 105 $i
                
                #python ms_network_no_m.py jitter/chess_large_vel6on105off${i} ${k} chess_large_vel6on105off${i} 65. 65. 6 $i
                #python ms_network_rest_modet_gmn_40.py jitter/chess_large_vel6on105off${i} ${k} chess_large_vel6on105off${i} 130. 130. 6 $i
                #python ms_network_inspect.py jitter/chess_size${q}_vel6on105off${i} ${k} chess_size${q}_vel6on105off${i} 131. 131. 105 $i ${m}.0
                
                #python video_evaluation.py jitter/chess_size${q}_vel6on105off${i} ${k} chess_size${q}_vel6on105off${i} 131. 131. 105 $i ${m}.0
                
                #python ms_network_no_m.py jitter/line_vel6on105off${i} ${k} line_vel6on105off${i} 65. 65. 6 $i
                #python ms_network_inspect.py jitter/line_vel6on105off${i} ${k} line_vel6on105off${i} 131. 131. 105 $i
                #python video_evaluation.py jitter/line_vel6on105off${i} ${k} line_vel6on105off${i} 131. 131. 105 $i
                
                python ms_network_midget_inspect.py murakami/on105off${i}/ ${k} on105off${i} 500. 140. 105 $i $m #481. 121.
                python video_evaluation.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i $m murakami
            done
        done 
    done
done


#bash network_poletti_rest_weight_modet_weight.sh
#bash together_poletti_rest.sh
#bash together_poletti_rest_modet.sh
#bash together_poletti_rest_weight_modet_weight.sh