#!/bin/bash

python PlaneProc.py -f  spectrums radardata fitting plotting -i ~/DATA/PerryPlane/MattsDataRed/ -c ~/DATA/PerryPlane/planeproc2red.ini  -r True
cp ~/DATA/PerryPlane/MattsDataRed/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages60secondsRed/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages60secondsRed/ -p ~/Dropbox/PerryPlane/fittedimages60secondsRed2/ -i 300
