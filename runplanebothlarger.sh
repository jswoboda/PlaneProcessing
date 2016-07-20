#!/bin/bash

python PlaneProc.py -f origdata spectrums radardata fitting plotting -m 4 -l 9 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_09 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_09/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_09/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Moving/Inputimages_09/

