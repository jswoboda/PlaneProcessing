#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PerryPlane/MattsDataInvfaster/ -c ~/DATA/PerryPlane/planeproc2red.ini -a ACF -k tikd -r True
cp -r ~/DATA/PerryPlane/MattsDataInvfaster/ACFInv/ ~/Dropbox/Inversion/Frameofref_faster/tikd
cp -r ~/DATA/PerryPlane/MattsDataInvfaster/FittedInv/ ~/Dropbox/Inversion/Frameofref_faster/tikd
cp -r ~/DATA/PerryPlane/MattsDataInvfaster/fittedimagestikd ~/Dropbox/Inversion/Frameofref_faster/tikd
