#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tik
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tikd
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tv