for k in 1 #2 3 4 5 6 7
do

    for v in 6 #0 3 6
    do
        for i in 0 15 30 45 60 75 90 105 
        do
            for q in 64 32 16 8 #4 2 96
            do
                echo chess_large off duty ${i}ms, velocity ${v}
                if [ ! -d ./video/img_input/jitter/chess_size${q}_vel${v}on105off${i}/${k} ]; then
                    mkdir -p ./video/img_input/jitter/chess_size${q}_vel${v}on105off${i}/${k}
                fi
                if [ ! -d ./data/jitter/chess_size${q}_vel${v}on105off${i}/${k} ]; then
                    mkdir -p ./data/jitter/chess_size${q}_vel${v}on105off${i}/${k}
                fi
                if [ ! -d ./data/jitter/chess_size${q}_vel${v}on105off${i}/network/${k} ]; then
                    mkdir -p ./data/jitter/chess_size${q}_vel${v}on105off${i}/network/${k};
                fi
                python image_creation_jitter.py chess_size${q}_vel${v}on105off${i}/${k} 105 $i $v 315 $q
        
        
#             echo rdp off duty ${i}ms, velocity ${v}
#             if [ ! -d ./video/img_input/jitter/rdp_vel${v}on105off${i}/${k} ]; then
#                 mkdir -p ./video/img_input/jitter/rdp_vel${v}on105off${i}/${k}
#             fi
#             if [ ! -d ./data/jitter/rdp_vel${v}on105off${i}/${k} ]; then
#                 mkdir -p ./data/jitter/rdp_vel${v}on105off${i}/${k}
#             fi
#             if [ ! -d ./data/jitter/rdp_vel${v}on105off${i}/network/${k} ]; then
#                 mkdir -p ./data/jitter/rdp_vel${v}on105off${i}/network/${k};
#             fi
#             python image_creation_jitter.py rdp_vel${v}on105off${i}/${k} 105 $i $v 315
#             
#             echo line off duty ${i}ms, velocity ${v}
#             if [ ! -d ./video/img_input/jitter/line_vel${v}on105off${i}/${k} ]; then
#                 mkdir -p ./video/img_input/jitter/line_vel${v}on105off${i}/${k}
#             fi
#             if [ ! -d ./data/jitter/line_vel${v}on105off${i}/${k} ]; then
#                 mkdir -p ./data/jitter/line_vel${v}on105off${i}/${k}
#             fi
#             if [ ! -d ./data/jitter/line_vel${v}on105off${i}/network/${k} ]; then
#                 mkdir -p ./data/jitter/line_vel${v}on105off${i}/network/${k};
#             fi
#             python image_creation_jitter.py line_vel${v}on105off${i}/${k} 105 $i $v 315
#         
            done
        done
    done
    

    echo creation done 
    
      for i in 0 15 30 45 60 75 90 105 
      do
          
          python ms_input.py jitter/chess_size8_vel6on105off${i}/${k} chess_size96_vel6on105off${i} 315 &
          python ms_input.py jitter/chess_size64_vel6on105off${i}/${k} chess_size64_vel6on105off${i} 315 &
          python ms_input.py jitter/chess_size32_vel6on105off${i}/${k} chess_size32_vel6on105off${i} 315 &
          python ms_input.py jitter/chess_size16_vel6on105off${i}/${k} chess_size16_vel6on105off${i} 315 
          wait
      done
    
#     python ms_input.py jitter/chess_vel0on105off0/${k} chess_vel0on105off0 315 &
#     python ms_input.py jitter/chess_vel0on105off30/${k} chess_vel0on105off30 315 &
#     python ms_input.py jitter/chess_vel0on105off45/${k} chess_vel0on105off45 315 &
#     python ms_input.py jitter/chess_vel0on105off60/${k} chess_vel0on105off60 315 
#     wait
#     
#     python ms_input.py jitter/chess_vel0on105off75/${k} chess_vel0on105off75 315 &
#     python ms_input.py jitter/chess_vel0on105off105/${k} chess_vel0on105off105 315 &
#     python ms_input.py jitter/chess_vel3on105off0/${k} chess_vel3on105off0 315 &
#     python ms_input.py jitter/chess_vel3on105off30/${k} chess_vel3on105off30 315 
#     wait
#     
#     python ms_input.py jitter/chess_vel3on105off45/${k} chess_vel3on105off45 315 &
#     python ms_input.py jitter/chess_vel3on105off60/${k} chess_vel3on105off60 315 &
#     python ms_input.py jitter/chess_vel3on105off75/${k} chess_vel3on105off75 315 &
#     python ms_input.py jitter/chess_vel3on105off105/${k} chess_vel3on105off105 315 
#     wait
#     
    #python ms_input.py jitter/chess_large_vel6on105off0/${k} chess_large_vel6on105off0 315 #&
    #python ms_input.py jitter/chess_large_vel6on105off0/${k} chess_large_vel6on105off0 315 #&
#     python ms_input.py jitter/chess_vel6on105off30/${k} chess_vel6on105off30 315 &
#     python ms_input.py jitter/chess_vel6on105off45/${k} chess_vel6on105off45 315 &
#     python ms_input.py jitter/chess_vel6on105off60/${k} chess_vel6on105off60 315 
#     wait
#     
#     python ms_input.py jitter/chess_vel6on105off75/${k} chess_vel6on105off75 315 &
#     python ms_input.py jitter/chess_vel6on105off105/${k} chess_vel6on105off105 315 &
#     python ms_input.py jitter/rdp_vel0on105off0/${k} rdp_vel0on105off0 315 &
#     python ms_input.py jitter/rdp_vel0on105off30/${k} rdp_vel0on105off30 315 
#     wait
#     
#     python ms_input.py jitter/rdp_vel0on105off45/${k} rdp_vel0on105off45 315 &
#     python ms_input.py jitter/rdp_vel0on105off60/${k} rdp_vel0on105off60 315 &
#     python ms_input.py jitter/rdp_vel0on105off75/${k} rdp_vel0on105off75 315 &
#     python ms_input.py jitter/rdp_vel0on105off105/${k} rdp_vel0on105off105 315 
#     wait
#     
#     python ms_input.py jitter/rdp_vel3on105off0/${k} rdp_vel3on105off0 315 &
#     python ms_input.py jitter/rdp_vel3on105off30/${k} rdp_vel3on105off30 315 &
#     python ms_input.py jitter/rdp_vel3on105off45/${k} rdp_vel3on105off45 315 &
#     python ms_input.py jitter/rdp_vel3on105off60/${k} rdp_vel3on105off60 315 
#     wait
#     
#     python ms_input.py jitter/rdp_vel3on105off75/${k} rdp_vel3on105off75 315 &
#     python ms_input.py jitter/rdp_vel3on105off105/${k} rdp_vel3on105off105 315 &
#     python ms_input.py jitter/rdp_vel6on105off0/${k} rdp_vel6on105off0 315 &
#     python ms_input.py jitter/rdp_vel6on105off30/${k} rdp_vel6on105off30 315 
#     wait
#     
#     python ms_input.py jitter/rdp_vel6on105off45/${k} rdp_vel6on105off45 315 &
#     python ms_input.py jitter/rdp_vel6on105off60/${k} rdp_vel6on105off60 315 &
#     python ms_input.py jitter/rdp_vel6on105off75/${k} rdp_vel6on105off75 315 &
#     python ms_input.py jitter/rdp_vel6on105off105/${k} rdp_vel6on105off105 315 
#     wait

    #python ms_input.py jitter/line_vel6on105off0/${k} line_vel6on105off0 315 &

 #   echo input done
         
    for v in 6 #0 3 6
    do
        for i in 0 15 30 45 60 75 90 105 
        do
            for q in 64 32 16 8 #4 2 96 
            do
                for m in 4 #2 4 6 8 
                do
                echo chess_large: $v $q
                #python graphics_chess.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $q
                python ms_network_no_m.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $m 200 #315
                #python ms_network_inspect.py jitter/chess_size${q}_vel6on105off${i} ${k} chess_size${q}_vel6on105off${i} 131. 131. 105 $i ${m}.0
                #python together_jitter_mdw.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $q $m
                #python ms_network_rest_no_gmn.py jitter/chess_large_vel${v}on105off${i} ${k} chess_large_vel${v}on105off${i} 130. 130. $v $i
                #python ms_network_rest_modet_gmn_40.py jitter/chess_large_vel${v}on105off${i} ${k} chess_large_vel${v}on105off${i} 130. 130. $v $i $m
                python ms_network_rest_modet_gmn_40.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 130. 130.  $v $i $m 
                python ms_network_inspect.py jitter/chess_size${q}_vel6on105off${i} ${k} chess_size${q}_vel6on105off${i} 131. 131. 105 $i ${m}.0
                python together_jitter_mdw.py jitter/chess_size${q}_vel${v}on105off${i} ${k} chess_size${q}_vel${v}on105off${i} 65. 65. $v $i $q $m
                #python video_evaluation.py jitter/chess_size${q}_vel6on105off${i} ${k} chess_size${q}_vel6on105off${i} 131. 131. 105 $i ${m}.0 jitter
                #python together_jitter.py jitter/chess_large_vel${v}on105off${i}/ ${k} chess_large_vel${v}on105off${i}
            
    #             echo rdp: $v $i
    #             python ms_network_no_m.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 65. 65. $v $i
    #             python ms_network_rest_no_gmn.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 130. 130. $v $i
    #             python ms_network_rest_modet_gmn_40.py jitter/rdp_vel${v}on105off${i} ${k} rdp_vel${v}on105off${i} 130. 130. $v $i
    #             python together_jitter.py jitter/rdp_vel${v}on105off${i}/ ${k} rdp_vel${v}on105off${i}

    #             echo line: $v $i
    #             python ms_network_no_m.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 65. 65. $v $i
    #             python ms_network_rest_no_gmn.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 130. 130. $v $i
    #             python ms_network_rest_modet_gmn_40.py jitter/line_vel${v}on105off${i} ${k} line_vel${v}on105off${i} 130. 130. $v $i
    #             python together_jitter.py jitter/line_vel${v}on105off${i}/ ${k} line_vel${v}on105off${i}
    #         
                done
            done
        done
    done
done
