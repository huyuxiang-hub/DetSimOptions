#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

if __name__ == "__main__":

    task = Sniper.TopTask("task")
    # task.asTop()
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
    tt_geom_svc = task.createSvc("TTGeomSvc")


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

    # = geometry service for PMT =
    import Geometry
    pmt_param_svc = task.createSvc("PMTParamSvc")
    tt_geom_svc = task.createSvc("TTGeomSvc")
   # tt_geom_svc.setLogLevel(6)

    # == DetSimOptions ==
    Sniper.loadDll("libDetSimOptions.so")
    Sniper.loadDll("libPhysiSim.so")
    Sniper.loadDll("libAnalysisCode.so")
    Sniper.loadDll("libCentralDetector.so")
    Sniper.loadDll("libTopTracker.so")
    Sniper.loadDll("libChimney.so")
    Sniper.loadDll("libCalibUnit.so")
    detsim0 = task.createSvc("DetSim0Svc")
    detsim0.property("AnaMgrList").set(["GenEvtInfoAnaMgr", "NormalAnaMgr"])
    #detsim0.property("PMTPosFile").set("Det1PMTPos_new.csv")
    detsim0.property("CDName").set("DetSim1")
    detsim0.property("TTName").set("")
    detsim0.property("VetoPMTPosMode").set("")
    print detsim0.objName()

    # == DetSimAlg ==
    Sniper.loadDll("libDetSimAlg.so")
    detsimalg = task.createAlg("DetSimAlg")
    detsimalg.property("DetFactory").set(detsim0.objName())
    detsimalg.property("VisMac").set("vis.mac")

    # === Central Detector ===
    cd = detsimalg.createTool("DetSim1Construction")
    cd.property("UseChimney").set(True)

    wp = detsimalg.createTool("WaterPoolConstruction")
    wp.property("enabledLatticedShell").set(True)

    # === latticed shell ===
    import math

    # ==== up xxx ====

    data_layer = [
                     ("up11", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2-6.0),
                     ("up10", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2),
                     ("up09", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0),
                     ("up08", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6),
                     ("up07", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8),
                     ("up06", 90.-8.4-8.4-8.4-8.4-8.2-8.0),
                     ("up05", 90.-8.4-8.4-8.4-8.4-8.2),
                     ("up04", 90.-8.4-8.4-8.4-8.4),
                     ("up03", 90.-8.4-8.4-8.4), 
                     ("up02", 90.-8.4-8.4), 
                     ("up01", 90.-8.4), 
                     ("equ", 90.),
                     ("bt01", 90.+8.4), 
                     ("bt02", 90.+8.4+8.4), 
                     ("bt03", 90.+8.4+8.4+8.4), 
                     ("bt04", 90.+8.4+8.4+8.4+8.4), # FIXME
                     ("bt05", 90.+8.4+8.4+8.4+8.4+8.2), 
                     ("bt06", 90.+8.4+8.4+8.4+8.4+8.2+8.0), 
                     ("bt07", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8), 
                     ("bt08", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6), 
                     ("bt09", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0), 
                     ("bt10", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2), 
                     ("bt11", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2+6.0), 
                     ]
    # build a dic
    dict_layer = dict(data_layer)
    print dict_layer

    data_up_down_30 = [
                     ("GLb3.up11", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2-6.0),
                     ("GLb4.up10", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0-6.2),
                     ("GLb3.up09", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6-7.0),
                     ("GLb2.up08", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8-7.6),
                     ("GLb2.up07", 90.-8.4-8.4-8.4-8.4-8.2-8.0-7.8),
                     ("GLb2.up06", 90.-8.4-8.4-8.4-8.4-8.2-8.0),
                     ("GLb1.up05", 90.-8.4-8.4-8.4-8.4-8.2),
                     ("GLb1.up04", 90.-8.4-8.4-8.4-8.4),
                     ("GLb1.up03", 90.-8.4-8.4-8.4), 
                     ("GLb1.up02", 90.-8.4-8.4), 
                     ("GLb1.up01", 90.-8.4), 
                     ("GLb2.equ", 90.),
                     ("GLb2.bt01", 90.+8.4), 
                     ("GLb1.bt02", 90.+8.4+8.4), 
                     ("GLb2.bt03", 90.+8.4+8.4+8.4), 
                     ("GLb2.bt04", 90.+8.4+8.4+8.4+8.4), # FIXME
                     ("GLb1.bt05", 90.+8.4+8.4+8.4+8.4+8.2), 
                     ("GLb1.bt06", 90.+8.4+8.4+8.4+8.4+8.2+8.0), 
                     ("GLb1.bt07", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8), 
                     ("GLb1.bt08", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6), 
                     ("GLb3.bt09", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0), 
                     ("GLb3.bt10", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2), 
                     ("GLb3.bt11", 90.+8.4+8.4+8.4+8.4+8.2+8.0+7.8+7.6+7.0+6.2+6.0), 
                     ]

    data_col_N = [
            ("GLw1.up10_up11", 6.0), 
            ("GLw1.up09_up10", 6.2), 
            ("GLw1.up08_up09", 7.0), 
            ("GLw1.up07_up08", 7.6), 
            ("GLw1.up06_up07", 7.8), 
            ("GLw1.up05_up06", 8.0), 
            ("GLw1.up04_up05", 8.2), 
            ("GLw1.up03_up04", 8.4), 
            ("GLw1.up02_up03", 8.4), 
            ("GLw1.up01_up02", 8.4), 
            ("GLw2.equ_up01",  8.4), 
            ("GLw2.equ_bt01",  8.4), 
            ("GLw3.bt01_bt02", 8.4), 
            ("GLw3.bt02_bt03", 8.4), 
            ("GLw2.bt03_bt04", 8.4), 
            ("GLw2.bt04_bt05", 8.2), 
            ("GLw1.bt05_bt06", 8.0), 
            ("GLw1.bt06_bt07", 7.8), 
            ("GLw1.bt07_bt08", 7.6), 
            ("GLw1.bt08_bt09", 7.0), 
            ("GLw1.bt09_bt10", 6.2), 
            ("GLw1.bt10_bt11", 6.0), 
            ]
    for label, theta in data_up_down_30:
        #l = 2*20050*math.sin(theta/180*math.pi)*math.sin(6./180*math.pi)
        GL = detsimalg.createTool("HBeamConstruction/%s"%label)
        #GL.property("L").set(l)

    for label, theta in data_col_N:
        #l = 2*20050*math.sin(0.5*theta/180*math.pi)
        GL = detsimalg.createTool("HBeamConstruction/%s"%label)
        #GL.property("L").set(l)

    # pillar
    data_pillar = [
            ("GZ1.A01_02", ()),
            ("GZ1.A02_03", ()),

            ("GZ1.A03_04", ()),
            ("GZ1.A04_05", ()),
            ("GZ1.A05_06", ()),
            ("GZ1.A06_07", ()),

            ("GZ1.B01_02", ()),
            ("GZ1.B02_03", ()),
            ("GZ1.B03_04", ()),
            ("GZ1.B04_05", ()),
            ("GZ1.B05_06", ()),
            ("GZ1.B06_07", ()),

            ("ZC2.A02_B02", ()),

            ("ZC2.A03_B03", ()),
            ("ZC2.A04_B04", ()),
            ("ZC2.A05_B05", ()),
            ("ZC2.A06_B06", ()),

            ("ZC2.A02_B03", ()),

            ("ZC2.A03_B04", ()),
            ("ZC2.A04_B05", ()),
            ("ZC2.A05_B06", ()),
            ("ZC2.A06_B07", ()),

            ("ZC2.B01_B01", ()),
            ("ZC2.B03_B03", ()),
            ("ZC2.B05_B05", ()),

            ("ZC2.A03_A03", ()),
            ("ZC2.A05_A05", ()),
            ]

    for label, theta in data_pillar:
        #l = 2*20050*math.sin(0.5*theta/180*math.pi)
        GL = detsimalg.createTool("HBeamConstruction/%s"%label)

    # ==== equator ====
    #enabled_list = ["GLb2.equ"]
    #if "GLb2.equ" in enabled_list:
    #    theta = 90.
    #    l = 2*20050*math.sin(theta/180*math.pi)*math.sin(6./180*math.pi)
    #    GL = detsimalg.createTool("HBeamConstruction/GLb2.equ")
    #    GL.property("L").set(l)

    # = begin run =
    task.show()
    task.run()
