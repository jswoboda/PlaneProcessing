#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Moving/Inputimages/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittederroronlyimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittederroronlyimages/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages/

python PlaneProc.py -f plotting -w y -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True

convert  ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittedimages/*.png ~/Dropbox/ISRErrorPapergifs/S1_2km_input.gif
convert  ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittedimages/*.png ~/Dropbox/ISRErrorPapergifs/S2_2km_output.gif 
