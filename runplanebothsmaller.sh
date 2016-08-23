#!/bin/bash


ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -m 6 -l 5 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_05 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_05/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Moving/Inputimages_05/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/fittederroronlyimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittederroronlyimages_05/

python PlaneProc.py -f  plotting -w y -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_05 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
convert ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/fittedimages/*.png ~/Dropbox/ISRErrorPapergifs/S4_10km_output.gif
