#!/bin/bash


ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red_statistics -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat_red_var.ini -r y

