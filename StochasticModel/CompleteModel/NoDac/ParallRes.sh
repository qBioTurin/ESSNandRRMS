#!/bin/bash
# 
#
# Usage: ./Resolution.sh  proc runs
#
#  1) solver 
#  2) type of Resolution
#  3) Initial Marking file


batch=$(( $2 / $1 ))


for r in $(eval echo "{1..$1}")
    do
    
    ri=$(( $batch*$r - $batch + 1 ))
    rf=$(($batch*$r ))
    
    nohup "./Resolution.sh" ModelRRMS.solver LSODA InitMark.txt $ri $rf &
    
done

