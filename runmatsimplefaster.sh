#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/SimpleDataFaster/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACFMat -k tikd -r False
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/SimpleDataFaster/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACFMat -k tik -r False
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/SimpleDataFaster/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACFMat -k tv -r False
#cp -r ~/DATA/PerryPlane/SimpleData/ACFInv/ ~/Dropbox/Inversion/Frameofrefsimple/tikd
#cp -r ~/DATA/PerryPlane/SimpleData/FittedInv/ ~/Dropbox/Inversion/Frameofrefsimple/tikd
#cp -r ~/DATA/PerryPlane/SimpleData/fittedimagestikd ~/Dropbox/Inversion/Frameofrefsimple/tikd
