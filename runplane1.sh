#!/bin/bash

python PlaneProc.py -f all -i ~/DATA/PlaneProcessing/exp_width_01 -c ~/DATA/PlaneProcessing/planeproc2.ini -r y
cp ~/DATA/PlaneProcessing/exp_width_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages/
