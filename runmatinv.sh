#!/bin/bash

ALLVAR=""
for var in "$@"
do
    ALLVAR="$ALLVAR $var"
done

if [ "$ALLVAR" = "" ]; then
    ALLVAR='plotting'
fi

python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tik -g .3
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/Tik/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/Fitted/*.h5 ~/Dropbox/Inversion/Tik/
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tikd -g .3
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/Tikd/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/Fitted/*.h5 ~/Dropbox/Inversion/Tikd/
python PlaneProcMat.py -f $ALLVAR -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -a ACF -k tv -g .3
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/ACFInv/*.h5 ~/Dropbox/Inversion/TV/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/Fitted/*.h5 ~/Dropbox/Inversion/TV/
