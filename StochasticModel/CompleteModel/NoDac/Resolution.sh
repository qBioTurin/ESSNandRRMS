#!/bin/bash
# 
#
# Usage: ./Resolution.sh  net.solver  type InitialMarking IDrunsInit IDrunsFinal
#
#  1) solver 
#  2) type of Resolution
#  3) Initial Marking file
#  4) runs
#------------------------------------------------------------------------
function float_gt() {
    perl -e "{if($1>$2){print 1} else {print 0}}"
}

###################################
##  Let define the injection times



    finaltime=8760

    InjEBVTime=(168 1608 3048 4488 5928)

    InjEBVTime+=( $finaltime ) 

    numberEBV=1000


    for r in $(eval echo "{$4..$5}")
        do
            nameSol="solution$r"
            nameRes="Resolution$r"
            name2inj="SecondInj.txt"
            nameInitM="InitialMarking$r.txt"
        
            echo 0 > "$name2inj"
            #  Into Solution.mtx should be stored the last time marking that can be used as initial marking for the next time window simulation.
            #
            start=`date +%s`
                
            ./$1 "$nameSol" -ftime $InjEBVTime[0]  -stime 24 -type $2 -taueps .05 -init $3 > log.txt

            cp "$nameSol.trace" "$nameRes.trace"

            END=$((${#InjEBVTime[@]}-1))

            for i in $(eval echo "{1..$END}")
            do     
            ## New final time of the simulation
            
                ii=$(($i-1))
                iii=$(($i+1))
                
            ## I need to save that the firs ebv inj is gone in order to start the mem activation
                if [ $i == 2 ] 
                    then echo 1 > $name2inj
                fi
                
                
                tempT=$(( ${InjEBVTime[$i]}-${InjEBVTime[$ii]} ))
            
            ## Read and modify the initial marking by adding the EBV cells
                
                InitialMarking=$( tail -n 1 "$nameSol.trace" )
                
                bkpIFS="$IFS"

                IFS=' ' read -r -a InitialMarking <<< $InitialMarking

                IFS="$bkpIFS"

                EBV=${InitialMarking[1]}

                InitialMarking[1]=$( echo $EBV + $numberEBV  | bc)
                
            ## I have to replace the negative values lower then 1e-9 otherwise the solver doesn't accept them!!

                a=1e-9
                b=0
                    
                for k in "${!InitialMarking[@]}"
                    do  
                    if [ $(float_gt $a ${InitialMarking[$k]}) == 1 ]
                    then InitialMarking[$k]=0
                    else if [  $(( $(float_gt  ${InitialMarking[$k]} $a) +  $(float_gt $b ${InitialMarking[$k]}) )) == 2 ] 
                            then  echo "ERROR: negative values!!!"
                        fi
                    fi
                    
                done         
                    
            ## I have to delete the time column

                unset InitialMarking[0]
                echo ${InitialMarking[@]} > $nameInitM
                
            ######################################
            ## Solution of the ODEs system or stoch. simulation
                
                ./$1 "$nameSol" -ftime $tempT  -stime 24  -type $2 -taueps .05 -init $nameInitM  > log.txt
                
            ## I have to add at the time column the initial time in order to have the right time points and delete the first row in the Resolution file beacasue when I merge the two files they have the same end and initial markings respectivaly.

                
                echo "$(tail -n +2 "$nameSol.trace")" > "$nameSol.trace"
                
                
                bkpIFS="$IFS"
                
                while  IFS=' ' read -r -a line
                do   
                    line[0]=$( echo ${line[0]} + ${InjEBVTime[$ii]} | bc)
                    echo ${line[@]} >> "temp$nameSol"
                    
                done < "$nameSol.trace" 
                
                cp "temp$nameSol" "$nameSol.trace"
                rm "temp$nameSol"
                
                IFS="$bkpIFS"
                
                temp="temp$r.txt"
                head -n -1 "$nameRes.trace" > $temp ;
                mv $temp "$nameRes.trace"
                
                cat "$nameSol.trace" >> ./"$nameRes.trace"
                
            done

            rm "$nameSol.trace"
            rm "$nameSol.mtx"
            rm $nameInitM
            rm $name2inj
            
            end=`date +%s`
            runtime=$(($end-$start))
            
            echo $runtime >> Runtime.txt
    done

