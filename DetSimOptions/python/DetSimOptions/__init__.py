#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao

import Sniper
# = geant4 related =
# == G4Svc ==
Sniper.loadDll("libG4Svc.so")
# == DetSimOptions ==
Sniper.loadDll("libDetSimOptions.so")
Sniper.loadDll("libAnalysisCode.so")
Sniper.loadDll("libPhysiSim.so")
Sniper.loadDll("libCentralDetector.so")
Sniper.loadDll("libTopTracker.so")
Sniper.loadDll("libChimney.so")
Sniper.loadDll("libCalibUnit.so")
# == DetSimAlg ==
Sniper.loadDll("libDetSimAlg.so")

# = helper =
# == load data from default directory ==
import os
def data_load(filename):
    # Load the data located in $DETSIMOPTIONSROOT/data
    detroot = os.environ["DETSIMOPTIONSROOT"]
    if not detroot:
        raise Exception("Can't get $DETSIMOPTIONSROOT")
    f = os.path.join(detroot, "data",  filename)
    if not os.path.exists(f):
        raise Exception("Can't get $DETSIMOPTIONSROOT/data/%s"%filename)

    return f

del Sniper
