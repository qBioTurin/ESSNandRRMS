#!/bin/bash
# 
#
# Usage: ./Resolution.sh  proc runs DAC
#
#  1) solver 
#  2) type of Resolution
#  3) if it is with DAC (1) or not (0)
#  4) solver

batch=$(( $2 / $1 ))


for r in $(eval echo "{1..$1}")
    do
    
    ri=$(( $batch*$r - $batch + 1 ))
    rf=$(($batch*$r ))
    
    nohup "./Resolution.sh" ModelRRMS.solver  $4 InitMark.txt $ri $rf $3 &
    
done

