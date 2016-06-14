#!/bin/bash

python PlaneProc.py -f all -i ~/DATA/PlaneProcessing/exp_width_stat_01 -c ~/DATA/PlaneProcessing/planeproc2_stat.ini -r y
cp ~/DATA/PlaneProcessing/exp_width_01_stat/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Stationary/fittedimages/
