#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

if __name__ == "__main__":

    task = Sniper.Task("task")
    task.asTop()
    task.setEvtMax(1)
    task.setLogLevel(0)
    # = I/O Related =
    import DataRegistritionSvc
    task.createSvc("DataRegistritionSvc")
    
    import RootIOSvc
    task.createSvc("RootOutputSvc/OutputSvc")
    ro = task.find("OutputSvc")
    ro.property("OutputStreams").set({"/Event/Sim": "sample_muon_without_op.root"})
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

    rootwriter.property("Output").set({"SIMEVT":"evt_muon_without_op.root"})

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun = gt.createTool("GtGunGenTool/gun")
    gun.property("particleNames").set(["mu-"])
    #gun.property("particleMomentums").set([215e3])
    gun.property("particleMomentums").set([100e3]) # 100 GeV
    #gun.property("DirectionMode").set("Fix")
    #gun.property("Directions").set([makeTV(0,0,0), makeTV(0,0,-1), makeTV(0,0,-1)])
    #gun.property("PositionMode").set("FixMany")
    gun.property("Positions").set([makeTV(0,0,0)])

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    import DetSimOptions
    from DetSimOptions.ConfAcrylic import ConfAcrylic
    acrylic_conf = ConfAcrylic(task)
    acrylic_conf.configure()
    acrylic_conf.set_gdml_output("geometry_acrylic.gdml")
    acrylic_conf.enable_chimney()
    acrylic_conf.set_tt_name("")

    acrylic_conf.add_anamgr("PhotonCollectAnaMgr")

    anamgr_gen = acrylic_conf.tool("GenEvtInfoAnaMgr")
    anamgr_gen.setLogLevel(6)
    anamgr_gen.property("EnableNtuple").set(True)
    anamgr_norm = acrylic_conf.tool("NormalAnaMgr")
    anamgr_norm.setLogLevel(6)
    anamgr_norm.property("EnableNtuple").set(True)
    # === PMT Mask ===
    #mask = acrylic_conf.mask()
    #print mask
    #mask.property("BufferMaterial").set("Water")
    # = begin run =
    task.show()
    task.run()
