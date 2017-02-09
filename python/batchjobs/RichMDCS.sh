#!/bin/bash

cd /afs/cern.ch/user/a/anandi/cmtuser/Brunel_v49r2p1/Rich/Rich1MDCS/python/

. SetupProject.sh Brunel v49r2p1

python runMDCSPeakFinder.py -R 171282 

echo 'Done'