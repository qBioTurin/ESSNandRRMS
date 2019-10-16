#!/bin/bash
# 
#
# Usage: ./Resolution.sh  proc runs DAC
#
#  1) solver 
#  2) type of Resolution
#  3) if it is with DAC or not


batch=$(( $2 / $1 ))


for r in $(eval echo "{1..$1}")
    do
    
    ri=$(( $batch*$r - $batch + 1 ))
    rf=$(($batch*$r ))
    
    nohup "./Resolution.sh" ModelRRMS.solver  SSA InitMark.txt $ri $rf $3 &
    
done

