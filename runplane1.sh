#!/bin/bash

python PlaneProc.py -f plotting -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_01 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -t 10 14 -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_01/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages/