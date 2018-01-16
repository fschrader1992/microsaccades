for k in 1 #2 3 4 5 6 7
do
    for v in 6 #0 3 6
    do
        for i in 0 #30 45 60 75 105
        do
#             echo chess: $v $i
#             python ms_network_no_m.py jitter/chess_vel${v}on105off${i} ${k} chess_vel${v}on105off${i} 65. 65. $v $i
#             python ms_network_rest_no_gmn.py jitter/chess_vel${v}on105off${i} ${k} chess_vel${v}on105off${i} 130. 130. $v $i
#             python ms_network_rest_modet_gmn_40.py jitter/chess_vel${v}on105off${i} ${k} chess_vel${v}on105off${i} 130. 130. $v $i
#             python together_jitter.py jitter/chess_vel${v}on105off${i}/ ${k} chess_vel${v}on105off${i}
#         
#             echo rdp: $v $i
#             python ms_network_no_m.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 65. 65. $v $i
#             python ms_network_rest_no_gmn.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 130. 130. $v $i
#             python ms_network_rest_modet_gmn_40.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 130. 130. $v $i
#             python together_jitter.py jitter/rdp_vel${v}on105off${i}/ ${k} rdp_vel${v}on105off${i}
#             
            echo line: $v $i
            python ms_network_no_m.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 65. 65. $v $i
            #python ms_network_rest_no_gmn.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 130. 130. $v $i
            #python ms_network_rest_modet_gmn_40.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 130. 130. $v $i
            #python together_jitter.py jitter/line_vel${v}on105off${i}/ ${k} line_vel${v}on105off${i}
        done
    done
done

