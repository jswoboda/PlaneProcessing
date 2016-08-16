#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -r y
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Stationary/Inputimages/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary/fittedimages/
