#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsData/ -c ~/DATA/PerryPlane/planeproc2.ini -r True
cp ~/DATA/PerryPlane/MattsData/Inputimages/*.png  ~/Dropbox/PerryPlane/Inputimages/
cp ~/DATA/PerryPlane/MattsData/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages60seconds/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages60seconds/ -p ~/Dropbox/PerryPlane/fittedimages60seconds2/ -i 60
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/Inputimages/ -p ~/Dropbox/PerryPlane/Inputimages2/ -i 15

python PlaneProc.py -f plotting -w y -i ~/DATA/PerryPlane/MattsData/ -c ~/DATA/PerryPlane/planeproc2.ini -r True
convert ~/DATA/PerryPlane/MattsData/Inputimages/*.png ~/Dropbox/ISRErrorPapergifs/S8_Output_FullParams.gif
