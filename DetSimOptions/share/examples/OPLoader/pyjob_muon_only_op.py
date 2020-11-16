#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

if __name__ == "__main__":

    Sniper.setLogLevel(6)
    task = Sniper.Task("task")
    task.asTop()
    task.setEvtMax(1)
    task.setLogLevel(6)
    # = I/O Related =
    import DataRegistritionSvc
    task.createSvc("DataRegistritionSvc")
    
    import RootIOSvc
    task.createSvc("RootOutputSvc/OutputSvc")
    ro = task.find("OutputSvc")
    ro.property("OutputStreams").set({"/Event/Sim": "sample_muon_only_op.root"})
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

    rootwriter.property("Output").set({"SIMEVT":"evt_muon_only_op.root"})

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    oploader = gt.createTool("GtOPLoaderTool")
    oploader.property("inputFile").set("evt_muon_without_op.root")
    oploader.property("ChunkSize").set(1024)
    oploader.property("ChunkIndex").set(100)


    gt.property("GenToolNames").set([oploader.objName()])

    # = geant4 related =
    import DetSimOptions
    from DetSimOptions.ConfAcrylic import ConfAcrylic
    acrylic_conf = ConfAcrylic(task)
    acrylic_conf.configure()
    acrylic_conf.set_gdml_output("geometry_acrylic.gdml")
    acrylic_conf.enable_chimney()
    acrylic_conf.set_tt_name("")

    #acrylic_conf.add_anamgr("PhotonCollectAnaMgr")

    acrylic_conf.detsimfactory().property("AnaMgrList").set(["DataModelWriter"])
    acrylic_conf.add_anamgr("NormalAnaMgr")

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
