#!/bin/bash

python PlaneProc.py -f radardata fitting plotting -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_state_01_red -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat_red.ini -r y
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary_Reduced/fittedimages/
