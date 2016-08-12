#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -m 4 -l 9 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_09 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_09/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/fittederrorimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittederrorimages_09/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Moving/Inputimages_09/


python PlaneProc.py -f  plotting -w y -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_09 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini

convert ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/fittedimages/*.png S6_18km_output.gif
