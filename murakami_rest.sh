for k in 11 #3 4 5 6 7 8 9 10 11
do
    for i in 0 15 30 45 60 75 90 105
    do
        python ms_network.py murakami/on105off${i} ${k} on105off${i} 481. 121. 105 $i
        python together_lr.py murakami/on105off${i} ${k} on105off${i} 481. 121. 105 $i 4
		
		#python ms_network_no_m.py murakami/on105off${i} ${k} on105off${i} 481. 121. 105 $i
        #python data_clean.py murakami/on105off${i}/ ${k} on105off${i}
        
		
		#python ms_network_rest.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        #python ms_network_rest_modet.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        #python ms_network_rest_weight.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        #python ms_network_rest_modet_weight.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        
		
		#python ms_network_rest_modet_gmn_40_murakami.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        
		
		#python together_murakami.py murakami/on105off${i}/ ${k} on105off${i}
        #python together_murakami_modet.py murakami/on105off${i}/ ${k} on105off${i}
        #python together_murakami_weight.py murakami/on105off${i}/ ${k} on105off${i}
        #python together_murakami_modet_weight.py murakami/on105off${i}/ ${k} on105off${i}
        
		
		#python together_murakami_modet_gmn_40.py murakami/on105off${i}/ ${k} on105off${i}
        
        
		#python ms_network_rest_no_gmn.py murakami/on105off${i}/ ${k} on105off${i} 481. 121. 105 $i
        #python together_murakami_no_gmn.py murakami/on105off${i}/ ${k} on105off${i}
    done
done