#!/bin/bash

python ms_input.py poletti2010/exp1/cond4_light/1_8_pi_arc/1 exp1_cond4_light_1_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/2_8_pi_arc/1 exp1_cond4_light_2_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/3_8_pi_arc/1 exp1_cond4_light_3_8_pi_arc_1 1000 
wait
date
echo second three files processed

python ms_input.py poletti2010/exp2/cond7/5_8_pi_arc/1 exp2_cond7_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/7_8_pi_arc/1 exp2_cond7_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/10_8_pi_arc/1 exp2_cond7_10_8_pi_arc_1 1000 
wait
date
echo third three files processed

python ms_input.py poletti2010/exp1/cond4_light/5_8_pi_arc/1 exp1_cond4_light_5_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/7_8_pi_arc/1 exp1_cond4_light_7_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/10_8_pi_arc/1 exp1_cond4_light_10_8_pi_arc_1 1000 
wait
date
echo fourth three files processed

python ms_input.py poletti2010/exp2/cond7/12_8_pi_arc/1 exp2_cond7_12_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/14_8_pi_arc/1 exp2_cond7_14_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp2/cond7/15_8_pi_arc/1 exp2_cond7_15_8_pi_arc_1 1000 
wait
date
echo fifth three files processed

python ms_input.py poletti2010/exp1/cond4_light/12_8_pi_arc/1 exp1_cond4_light_12_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/14_8_pi_arc/1 exp1_cond4_light_14_8_pi_arc_1 1000 &
python ms_input.py poletti2010/exp1/cond4_light/15_8_pi_arc/1 exp1_cond4_light_15_8_pi_arc_1 1000 
wait
date
echo sixth three files processed