#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao


import Sniper

def get_parser():

    import argparse

    parser = argparse.ArgumentParser(description='Run Prototype Simulation.')

    parser.add_argument("--evtmax", type=int, default=100, help='events to be processed')
    parser.add_argument("--seed", type=int, default=42, help='seed')
    parser.add_argument("--vis", action="store_true", default=False)
    parser.add_argument("--output", default="sample_detsim.root", help="output file name")
    parser.add_argument("--user-output", default="sample_detsim_user.root", help="output file name")
    parser.add_argument("--mac", default="run.mac", help="mac file")
    parser.add_argument("--particle", default="gamma", help="particle name")
    parser.add_argument("--momentum", default=1., type=float, help="momentum size (MeV)")

    return parser

if __name__ == "__main__":

    parser = get_parser()
    args = parser.parse_args()

#############################################################################
    task = Sniper.Task("task")
    task.setEvtMax(args.evtmax)
    #task.setLogLevel(0)
    # = I/O Related =
    import DataRegistritionSvc
    task.createSvc("DataRegistritionSvc")
    
    import RootIOSvc
    task.createSvc("RootOutputSvc/OutputSvc")
    ro = task.find("OutputSvc")
    ro.property("OutputStreams").set({"/Event/Sim": args.output})
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

    rootwriter.property("Output").set({"SIMEVT":args.user_output})

    # = generator related =
    import GenTools
    from GenTools import makeTV
    gt = task.createAlg("GenTools")

    gun = gt.createTool("GtGunGenTool/gun")
    gun.property("particleNames").set([args.particle])
    gun.property("particleMomentums").set([args.momentum])
    #gun.property("DirectionMode").set("Fix")
    #gun.property("Directions").set([makeTV(0,0,0), makeTV(0,0,-1), makeTV(0,0,-1)])
    #gun.property("PositionMode").set("FixMany")
    gun.property("Positions").set([makeTV(0,0,0)])

    gt.property("GenToolNames").set([gun.objName()])

    # = geant4 related =
    import DetSimOptions
    g4svc = task.createSvc("G4Svc")

        # == DetSimOptions ==
    factory = task.createSvc("DetSim0Svc")
    factory.property("CDName").set("Prototype")
    factory.property("PMTMother").set("lBuffer")
    factory.property("TTName").set("")
    factory.property("AnaMgrList").set(["DataModelWriter",
                                        "GenEvtInfoAnaMgr",
                                        "NormalAnaMgr", 
                                        "DepositEnergyAnaMgr",
                                        "InteresingProcessAnaMgr"])
    detsimalg = task.createAlg("DetSimAlg")
    detsimalg.property("DetFactory").set(factory.objName())
    detsimalg.property("RunMac").set("run.mac")
    if args.vis:
        detsimalg.property("VisMac").set("vis.mac")
    #detfactory.property("3inchPMTPosFile").set("")
    # === Prototype Construction ===
    prototype = detsimalg.createTool("PrototypeConstruction")
    prototype.property("UseEquatorCircle").set(True)
    prototype.property("UseChimney").set(True)
    # === PMT Mask ===
    #mask = acrylic_conf.mask()
    #print mask
    #mask.property("BufferMaterial").set("Water")
    pmtsdmgr = detsimalg.createTool("PMTSDMgr")
    pmtsdmgr.property("PMTSD").set("dywSD_PMT_v2")
    op = detsimalg.createTool("DsPhysConsOptical")
    op.property("UsePMTOpticalModel").set(False)
    # = begin run =
    task.show()
    task.run()
