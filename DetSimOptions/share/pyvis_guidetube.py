#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

if __name__ == "__main__":

    task = Sniper.Task("task")
    task.setEvtMax(10)
    task.setLogLevel(0)
    # = Data Buffer =
    import BufferMemMgr
    bufMgr = task.createSvc("BufferMemMgr")

    # = random svc =
    import RandomSvc
    task.property("svcs").append("RandomSvc")
    rndm = task.find("RandomSvc")
    rndm.property("Seed").set(42)

    # = root writer =
    import RootWriter
    rootwriter = task.createSvc("RootWriter")

    rootwriter.property("Output").set({"SIMEVT":"evt.root"})

    # = geometry service for PMT =
    import Geometry
    pmt_param_svc = task.createSvc("PMTParamSvc")

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun = gt.createTool("GtGunGenTool/gun")
    gun.property("particleNames").set(["gamma", "gamma", "e+"])
    gun.property("particleMomentums").set([1.173, 1.333, 1.0])
    gun.property("DirectionMode").set("Fix")
    gun.property("Directions").set([makeTV(0,0,-1), makeTV(0,0,-1), makeTV(0,0,-1)])
    gun.property("PositionMode").set("FixMany")
    gun.property("Positions").set([makeTV(0,0,15500), makeTV(0,0,15500), makeTV(0,0,-155)])

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    # == G4Svc ==
    Sniper.loadDll("libG4Svc.so")
    g4svc = task.createSvc("G4Svc")
    print g4svc.objName()

    # == DetSimOptions ==
    Sniper.loadDll("libDetSimOptions.so")
    Sniper.loadDll("libAnalysisCode.so")
    Sniper.loadDll("libCentralDetector.so")
    Sniper.loadDll("libTopTracker.so")
    Sniper.loadDll("libChimney.so")
    Sniper.loadDll("libCalibUnit.so")
    detsim0 = task.createSvc("DetSim0Svc")
    detsim0.property("AnaMgrList").set(["GenEvtInfoAnaMgr", "NormalAnaMgr"])
    detsim0.property("PMTPosFile").set("PMTPos_Acrylic_with_chimney.csv")
    detsim0.property("CDName").set("DetSim1")
    detsim0.property("TTName").set("")
    detsim0.property("VetoPMTPosMode").set("")
    print detsim0.objName()

    # == DetSimAlg ==
    Sniper.loadDll("libDetSimAlg.so")
    detsimalg = task.createAlg("DetSimAlg")
    detsimalg.property("DetFactory").set(detsim0.objName())
    detsimalg.property("VisMac").set("vis_guidetube.mac")

    # === Central Detector ===
    cd = detsimalg.createTool("DetSim1Construction")
    cd.property("UseChimney").set(True)
    cd.property("CheckOverlap").set(True)

    lowerchim = detsimalg.createTool("LowerChimney")
    lowerchim.property("UseShutter").set(False)

    # === Guide Tube ===
    detsim0.property("CalibUnitEnable").set(True)
    detsim0.property("CalibUnitName").set("Calib_GuideTube_V1") # this is default one
    detsim0.property("CalibUnitExtras").set([#"Calib_GuideTube_V1_0", 
                                             #"Calib_GuideTube_V1_1"
                                             ]) # Enable more calibration units here

    guide_tube_v1_0 = detsimalg.createTool("Calib_GuideTube_Construction/Calib_GuideTube_Construction_V1_0")
    guide_tube_v1_0.property("Option").set("V1_0")

    guide_tube_v1_1 = detsimalg.createTool("Calib_GuideTube_Construction/Calib_GuideTube_Construction_V1_1")
    guide_tube_v1_1.property("Option").set("V1_1")

    import math

    guide_tube_place_v1_0 = detsimalg.createTool("Calib_GuideTube_Placement/Calib_GuideTube_Placement_V1_0")
    guide_tube_place_v1_0.property("Phi").set(-(90.+2.6)*math.pi/180.)

    guide_tube_place_v1_1 = detsimalg.createTool("Calib_GuideTube_Placement/Calib_GuideTube_Placement_V1_1")
    guide_tube_place_v1_1.property("Phi").set((303.4-180)*math.pi/180.)
    # = begin run =
    task.show()
    task.run()
