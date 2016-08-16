#!/bin/bash


ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat_red.ini -r y
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary_Reduced/fittedimages/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Stationary_Reduced/Inputimages/

