#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

def get_parser():

    import argparse

    parser = argparse.ArgumentParser(description='Run Atmospheric Simulation.')
    parser.add_argument("--evtmax", type=int, default=100, help='events to be processed')
    parser.add_argument("--seed", type=int, default=42, help='seed')
    parser.add_argument("--output", default="evt.root", help="output file name")
    parser.add_argument("--mac", default="run.mac", help="mac file")
    parser.add_argument("--gst", default="/afs/ihep.ac.cn/users/g/guowl/genie/105.root", help="gst file")
    parser.add_argument("--gstindex", type=int, default=0, help='gst start index')

    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    task = Sniper.TopTask("task")
    
    task.setEvtMax(args.evtmax)
    #task.setEvtMax(25)
    #task.setLogLevel(0)
    task.setLogLevel(5)
    # = Data Buffer =
    import BufferMemMgr
    bufMgr = task.createSvc("BufferMemMgr")

    # = random svc =
    import RandomSvc
    task.property("svcs").append("RandomSvc")
    rndm = task.find("RandomSvc")
    rndm.property("Seed").set(args.seed)

    # = root writer =
    import RootWriter
    rootwriter = task.createSvc("RootWriter")

    rootwriter.property("Output").set({"SIMEVT":args.output})

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun_gst = gt.createTool("GtGstTool")
    gun_gst.property("inputGstFile").set(args.gst)
    gun_gst.property("GstStartIndex").set(args.gstindex)

    gun_pos = gt.createTool("GtPositionerTool")
    gun_pos.property("GenInVolume").set("pTarget")
    gun_pos.property("Material").set("LS")

    gt.property("GenToolNames").set([gun_gst.objName(), gun_pos.objName()])

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
    detsim0 = task.createSvc("DetSim0Svc")
    detsim0.property("AnaMgrList").set(["GenEvtInfoAnaMgr", "NormalAnaMgr",
                                        "DepositEnergyAnaMgr", 
                                        "InteresingProcessAnaMgr"])
    detsim0.property("PhysicsList").set("QGSP_BERT_HP")
    detsim0.property("CDName").set("DetSim1")
    detsim0.property("TTName").set("")
    print detsim0.objName()

    # == DetSimAlg ==
    Sniper.loadDll("libDetSimAlg.so")
    detsimalg = task.createAlg("DetSimAlg")
    detsimalg.property("DetFactory").set(detsim0.objName())
    detsimalg.property("RunMac").set(args.mac)
    detsimalg.property("RunCmds").set([
            "/run/initialize",
            "/process/inactivate Scintillation",
            "/process/inactivate Cerenkov",
        ])
    # === configure GenEvtInfoAnaMgr ===
    anamgr_genevtinfo = detsimalg.createTool("GenEvtInfoAnaMgr")
    anamgr_genevtinfo.property("EnableNtuple").set(True)
    # === configure NormalAnaMgr ===
    anamgr_normal = detsimalg.createTool("NormalAnaMgr")
    anamgr_normal.property("EnableNtuple").set(True)
    # === configure DepositEnergyAnaMgr ===
    anamgr_edep = detsimalg.createTool("DepositEnergyAnaMgr")
    anamgr_edep.property("EnableNtuple").set(True)

    # = begin run =
    task.show()
    task.run()
