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
    parser.add_argument("--gen-x", default=0.0, type=float)
    parser.add_argument("--gen-y", default=0.0, type=float)
    parser.add_argument("--gen-z", default=0.0, type=float)
    parser.add_argument("--dir-x", default=0.0, type=float)
    parser.add_argument("--dir-y", default=0.0, type=float)
    parser.add_argument("--dir-z", default=0.0, type=float)
    parser.add_argument("--split-mode", default="PrimaryTrack", 
                                   choices=["PrimaryTrack",
                                            "EveryTrack",
                                            "Time"
                                            ])
    parser.add_argument("--split-time", default=3000., type=float)
    parser.add_argument("--voxel-merge-flag", action="store_true", dest="voxel_merge_flag")
    parser.add_argument("--voxel-no-merge-flag", action="store_false", dest="voxel_merge_flag")
    parser.set_defaults(voxel_merge_flag=False)
    parser.add_argument("--voxel-merge-twin", default=1, type=float)
    parser.add_argument("--voxel-fill-ntuple", action="store_true", dest="voxel_fill_ntuple")
    parser.add_argument("--voxel-no-fill-ntuple", action="store_true", dest="voxel_fill_ntuple")
    parser.set_defaults(voxel_fill_ntuple=False)

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
    #gt.setLogLevel(2)

    gun = gt.createTool("GtGunGenTool/gun")
    #gun.property("particleNames").set(["mu-"])
    #gun.property("particleMomentums").set([215e3])
    #gun.property("particleMomentums").set([100e3]) # 100 GeV
    if args.dir_x != 0.0 or args.dir_y != 0.0 or args.dir_z != 0:
        gun.property("DirectionMode").set("Fix")
        gun.property("Directions").set([makeTV(args.dir_x,args.dir_y,args.dir_z)])
    gun.property("PositionMode").set("FixMany")
    from GenTools import makeTV
    gun.property("Positions").set([makeTV(args.gen_x, args.gen_y, args.gen_z)])

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

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    import DetSimOptions
    from DetSimOptions.ConfAcrylic import ConfAcrylic
    acrylic_conf = ConfAcrylic(task)
    acrylic_conf.configure()
    acrylic_conf.set_tt_name("")
    acrylic_conf.enable_PMTSD_v2()
    acrylic_conf.disable_pmts_and_struts_in_cd()

    # voxel method
    # check the data path
    import os
    if os.environ.get("DETSIMOPTIONSROOT", None) is None:
        print "Missing DETSIMOPTIONSROOT"
        sys.exit(-1)
    dp = lambda f: os.path.join(os.environ.get("DETSIMOPTIONSROOT"),
                                "share", "examples", "voxelmuon", f)
    acrylic_conf.add_anamgr("MuonFastSimVoxel")
    mfsv = acrylic_conf.tool("MuonFastSimVoxel")
    mfsv.property("GeomFile").set(dp("geom-geom-20pmt.root"))
    mfsv.property("NPEFile").set(dp("npehist3d_single.root"))
    mfsv.property("HitTimeMean").set(dp("hist3d.root"))
    mfsv.property("HitTimeRes").set(dp("dist_tres_single.root"))
    mfsv.property("MergeFlag").set(args.voxel_merge_flag)
    mfsv.property("MergeTimeWindow").set(args.voxel_merge_twin)
    mfsv.property("EnableNtuple").set(args.voxel_fill_ntuple)

    acrylic_conf.add_anamgr("PostponeTrackAnaMgr")

    pta = acrylic_conf.tool("PostponeTrackAnaMgr")
    #pta.setLogLevel(2)
    pta.property("SplitMode").set(args.split_mode)
    pta.property("TimeCut").set(args.split_time)

    op_process = acrylic_conf.optical_process()
    op_process.property("UseScintillation").set(False)
    op_process.property("UseCerenkov").set(False)

    # = begin run =
    task.show()
    task.run()
