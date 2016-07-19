#!/bin/bash

python PlaneProc.py -f plotting -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2_stat.ini -r y
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_stat_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary/fittedimages/
