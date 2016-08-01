#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProc.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsDataRed/ -c ~/DATA/PerryPlane/planeproc2red.ini  -r True
cp ~/DATA/PerryPlane/MattsDataRed/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages60secondsRed/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages60secondsRed/ -p ~/Dropbox/PerryPlane/fittedimages60secondsRed2/ -i 420
