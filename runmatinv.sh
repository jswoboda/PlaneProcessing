#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tik -g 0.1
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/Tik/ACFInv
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/FittedInv/*.h5 ~/Dropbox/Inversion/Tik/FittedInv
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tikd -g 2
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/Tikd/ACFInv
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/FittedInv/*.h5 ~/Dropbox/Inversion/Tikd/FittedInv
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tv -g 0.02
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/TV/ACFInv
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/FittedInv/*.h5 ~/Dropbox/Inversion/TV/FittedInv
