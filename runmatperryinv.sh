#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi


python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsDataInv/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACF -k tik 
cp -r ~/DATA/PerryPlane/MattsDataInv/ACFInv/ ~/Dropbox/Inversion/Frameofref/tik
cp -r ~/DATA/PerryPlane/MattsDataInv/FittedInv/ ~/Dropbox/Inversion/Frameofref/tik
cp -r ~/DATA/PerryPlane/MattsDataInv/fittedimagestik ~/Dropbox/Inversion/Frameofref/tik

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsDataInv/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACF -k tikd 
cp -r ~/DATA/PerryPlane/MattsDataInv/ACFInv/ ~/Dropbox/Inversion/Frameofref/tikd
cp -r ~/DATA/PerryPlane/MattsDataInv/FittedInv/ ~/Dropbox/Inversion/Frameofref/tikd
cp -r ~/DATA/PerryPlane/MattsDataInv/fittedimagestikd ~/Dropbox/Inversion/Frameofref/tikd

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsDataInv/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACF -k tv 
cp -r ~/DATA/PerryPlane/MattsDataInv/ACFInv/ ~/Dropbox/Inversion/Frameofref/tv
cp -r ~/DATA/PerryPlane/MattsDataInv/FittedInv/ ~/Dropbox/Inversion/Frameofref/tv
cp -r ~/DATA/PerryPlane/MattsDataInv/fittedimagestv ~/Dropbox/Inversion/Frameofref/tv
