#!/bin/bash

python PlaneProc.py -f radardata fitting plotting -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat_red.ini -r y
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01_red/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary_Reduced/fittedimages/
