#!/bin/bash

python PlaneProc.py -f plotting -i ~/DATA/PerryPlane/MattsData/ -c ~/DATA/PerryPlane/planeproc2.ini -r y
cp ~/DATA/PerryPlane/MattsData/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages60seconds/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages60seconds/ -p ~/Dropbox/PerryPlane/fittedimages60seconds2/ -i 60
