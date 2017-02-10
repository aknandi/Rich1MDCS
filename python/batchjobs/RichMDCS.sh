#!/bin/bash

username=$USER

source /afs/cern.ch/lhcb/software/releases/LBSCRIPTS/dev/InstallArea/scripts/LbLogin.sh

cd /afs/cern.ch/user/${username:0:1}/$username/cmtuser/BrunelDev_v52r0/
 
./run python Rich/Rich1MDCS/python/runMDCSPeakFinder.py -R 171282 

echo 'Done'