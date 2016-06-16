#!/bin/bash

python PlaneProc.py -f radardata fitting plotting -i ~/DATA/PerryPlane/MattsData240sec/ -c ~/DATA/PerryPlane/planeproc240sec.ini -t 2 -r True
cp ~/DATA/PerryPlane/MattsData240sec/Inputimages/*.png  ~/Dropbox/PerryPlane/Inputimages/
cp ~/DATA/PerryPlane/MattsData240sec/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages240seconds/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages240seconds/ -p ~/Dropbox/PerryPlane/fittedimages240seconds2/ -i 240
