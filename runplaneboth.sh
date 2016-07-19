#!/bin/bash

python PlaneProc.py -f origdata spectrums radardata fitting plotting -m 6 -l 3 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_03 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_03/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_03/

python PlaneProc.py -f origdata spectrums radardata fitting plotting -m 3 -l 5 -i ~/DATA/PlaneProcessing/Phased_Array/exp_width_06 -c ~/DATA/PlaneProcessing/Phased_Array/planeproc2.ini -r True
cp ~/DATA/PlaneProcessing/Phased_Array/exp_width_06/fittedimages/*.png  ~/Dropbox/PlaneProcessing/Moving/fittedimages_06/
