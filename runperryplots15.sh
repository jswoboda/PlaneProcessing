#!/bin/bash

python PlaneProc.py -f radardata fitting plotting -i ~/DATA/PerryPlane/MattsData15sec/ -c ~/DATA/PerryPlane/planeproc15sec.ini -r n
cp ~/DATA/PerryPlane/MattsData15sec/fittedimages/*.png  ~/Dropbox/PerryPlane/fittedimages15seconds/
python PlaneProcUtil.py -f ~/Dropbox/PerryPlane/fittedimages15seconds/ -p ~/Dropbox/PerryPlane/fittedimages15seconds2/ -i 15
