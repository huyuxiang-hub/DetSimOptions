#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

def get_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Run Detector Simulation.')
    parser.add_argument("--evtmax", type=int, default=1, help='events to be processed')
    parser.add_argument("--types", default="e-", 
                                   choices=["e-",
                                            "mu-",
                                            ])
    parser.add_argument("--split-mode", default="PrimaryTrack", 
                                   choices=["PrimaryTrack",
                                            "EveryTrack",
                                            "Time"
                                            ])
    parser.add_argument("--split-time", default=3000., type=float)

    return parser

if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()

    task = Sniper.Task("task")
    task.asTop()
    task.setEvtMax(args.evtmax)
    task.setLogLevel(3)
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
    gt.setLogLevel(2)

    gun = gt.createTool("GtGunGenTool/gun")
    #gun.property("particleNames").set(["mu-"])
    #gun.property("particleMomentums").set([215e3])
    #gun.property("particleMomentums").set([100e3]) # 100 GeV
    #gun.property("DirectionMode").set("Fix")
    #gun.property("Directions").set([makeTV(0,0,0), makeTV(0,0,-1), makeTV(0,0,-1)])
    #gun.property("PositionMode").set("FixMany")

    # muon
    if args.types == "mu-":
        gun.property("particleNames").set(["mu-"])
        gun.property("particleMomentums").set([100e3]) # 100 GeV
    elif args.types == "e-":
    # electron
        gun.property("particleNames").set(["e-"])
        gun.property("particleMomentums").set([1.])
    else:
        print "Unknown types"
        import sys
        sys.exit(-1)

    gun.property("Positions").set([makeTV(0,0,0)])

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    import DetSimOptions
    from DetSimOptions.ConfAcrylic import ConfAcrylic
    acrylic_conf = ConfAcrylic(task)
    acrylic_conf.configure()
    acrylic_conf.set_tt_name("")
    acrylic_conf.enable_PMTSD_v2()
    acrylic_conf.disable_pmts_and_struts_in_cd()

    acrylic_conf.add_anamgr("PostponeTrackAnaMgr")

    pta = acrylic_conf.tool("PostponeTrackAnaMgr")
    pta.setLogLevel(2)
    pta.property("SplitMode").set(args.split_mode)
    pta.property("TimeCut").set(args.split_time)

    op_process = acrylic_conf.optical_process()
    op_process.property("UseScintillation").set(False)
    op_process.property("UseCerenkov").set(False)

    # = begin run =
    task.show()
    task.run()
