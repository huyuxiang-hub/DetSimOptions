#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

if __name__ == "__main__":

    task = Sniper.Task("task")
    task.asTop()
    task.setEvtMax(3)
    task.setLogLevel(0)
    # = I/O Related =
    import DataRegistritionSvc
    task.createSvc("DataRegistritionSvc")
    
    import RootIOSvc
    task.createSvc("RootOutputSvc/OutputSvc")
    ro = task.find("OutputSvc")
    ro.property("OutputStreams").set({"/Event/Sim": "sample.root"})
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

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun = gt.createTool("GtGunGenTool/gun")
    gun.property("particleNames").set(["gamma"])
    gun.property("particleMomentums").set([1.0])
    #gun.property("DirectionMode").set("Fix")
    #gun.property("Directions").set([makeTV(0,0,0), makeTV(0,0,-1), makeTV(0,0,-1)])
    #gun.property("PositionMode").set("FixMany")
    gun.property("Positions").set([makeTV(0,0,0)])

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    import DetSimOptions
    from DetSimOptions.ConfBalloon import ConfBalloon
    balloon_conf = ConfBalloon(task)
    balloon_conf.configure()
    balloon_conf.set_gdml_output("geometry_balloon.gdml")
    
    detfactory = balloon_conf.detsimfactory()
    #detfactory.property("3inchPMTPosFile").set("")
    #detfactory.property("PMTName").set("PMTMask")
    # === PMT Mask ===
    #mask = balloon_conf.mask()
    #print mask
    #mask.property("BufferMaterial").set("LAB")
    # === Ballon CD ===
    ballooncd = balloon_conf.ballooncd()
#   ballooncd.property("BalloonMaterial").set("PA")
#   ballooncd.property("BalloonThickness").set(0.1)

    # = begin run =
    task.show()
    task.run()
