#!/bin/bash
for k in 1 2 3 4 5 6 7
do
    for i in 0 15 30 45 60 75 90 105
    do
    
        python together_murakami.py murakami/on105off${i}/ ${k} on105off${i}
    
    done
done