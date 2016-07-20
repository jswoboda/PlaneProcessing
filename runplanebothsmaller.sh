#!/bin/bash

python PlaneProc.py -f origdata spectrums radardata fitting plotting -m 6 -l 5 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_05 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_05/
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_05/Inputimages/*.png  ~/Dropbox/PlaneProcessing/Moving/Inputimages_05/
