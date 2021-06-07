#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao

import os.path
import Sniper
import DetSimOptions

class ConfAcrylic(object):

    def __init__(self, task):
        self._task = task

        self._g4svc = None
        self._factory = None
        self._detsimalg = None
        self._mask = None

        self._geom = None
        self._gdml = None
        self._dae = None

        self._pmtsd_mgr = None
        self._optical_process = None
        self._em_process = None

    def configure(self):
        # == G4Svc ==
        self._g4svc = self._task.createSvc("G4Svc")

        # == DetSimOptions ==
        self._factory = self._task.createSvc("DetSim0Svc")
        self._factory.property("AnaMgrList").set(["DataModelWriter",
                                                  "GenEvtInfoAnaMgr",
                                                  "NormalAnaMgr", 
                                                  "DepositEnergyAnaMgr",
                                                  "InteresingProcessAnaMgr"])
        self._factory.property("CDName").set("DetSim1")
        self._factory.property("PMTPosFile").set(DetSimOptions.data_load("PMTPos_Acrylic_with_chimney.csv"))
        self._factory.property("3inchPMTPosFile").set(DetSimOptions.data_load("3inch_pos.csv"))
        self._factory.property("StrutPosFile").set(DetSimOptions.data_load("Strut_Acrylic.csv"))
        self._factory.property("Strut2PosFile").set(DetSimOptions.data_load("StrutBar2_Acrylic.csv"))
        self._factory.property("FastenerPosFile").set(DetSimOptions.data_load("Strut_Anchor_Acrylic.csv"))
        self._factory.property("XJanchorPosFile").set(DetSimOptions.data_load("XJanchor.csv"))
        self._factory.property("SJCLSanchorPosFile").set(DetSimOptions.data_load("SJCLSanchor.csv"))
        self._factory.property("SJReceiverPosFile").set(DetSimOptions.data_load("SJReceiverPos.csv"))
        self._factory.property("SJFixturePosFile").set(DetSimOptions.data_load("SJFixturePos.csv"))

        self._factory.property("ElecSimPMTQEScale").set(1.0) # or 0.9

        # == DetSimAlg ==
        self._detsimalg = self._task.createAlg("DetSimAlg")
        self._detsimalg.property("DetFactory").set(self._factory.objName())
        self._detsimalg.property("RunMac").set("run.mac")
        
    def detsimalg(self):
        if not self._detsimalg:
            self.configure()

        return self._detsimalg

    def detsimfactory(self):
        if not self._factory:
            self.configure()

        return self._factory

    def tool(self, toolname):
        if not self._detsimalg:
            self.configure()
        if self._detsimalg.findTool(toolname):
            return self._detsimalg.findTool(toolname)

        return self._detsimalg.createTool(toolname)

    def set_gdml_output(self, filename):
        # If already setup, just return
        if self._gdml and self._gdml != filename:
            raise Exception("You have already set the gdml output.")
        self._gdml = filename
        # flag to remove the temporary gdml, and don't generate the  
        # user wanted gdml file. 
        flag_should_remove = False 
        # if already exists, skip this
        if os.path.exists(self._gdml):
            print("The %s already exists, skip generate" %self._gdml)
            flag_should_remove = True
        # configure the GeoAnaMgr 
        geo = self.geom()
        geo.property("GdmlEnable").set(True)
        geo.property("GdmlOutput").set(self._gdml)
        geo.property("Gdml2RootEnable").set(True)
        geo.property("GdmlDeleteFile").set(flag_should_remove)
    def set_dae_output(self, filename):
        # If already setup, just return
        if self._dae and self._dae != filename:
            raise Exception("You have already set the dae output.")
        self._dae = filename
        # flag to remove the temporary dae, and don't generate the  
        # user wanted dae file. 
        flag_should_remove = False 
        # if already exists, skip this
        if os.path.exists(self._dae):
            print("The %s already exists, skip generate" %self._dae)
            return
        # configure the GeoAnaMgr 
        geo = self.geom()
        geo.property("DAEEnable").set(True)
        geo.property("DAEOutput").set(self._dae)

    def add_anamgr(self, toolname):
        self._factory.property("AnaMgrList").append(toolname)

    def set_muon_output(self):
        self._factory.property("AnaMgrList").append("MuProcessAnaMgr")
    def set_muonIso_output(self):
        self._factory.property("AnaMgrList").append("MuIsoProcessAnaMgr")
    def set_muonFastn_output(self):
        self._factory.property("AnaMgrList").append("MuFastnProcessAnaMgr")

    def set_tt_edep_output(self):
        self._factory.property("AnaMgrList").append("DepositEnergyTTAnaMgr")
        self.tt_anamgr = self.tool("DepositEnergyTTAnaMgr")
        self.tt_anamgr.property("EnableNtuple").set(True)

    def disable_pmts_and_struts_in_cd(self):
        self._factory.property("PMTPosFile").set("")
        self._factory.property("StrutPosFile").set("")
        self._factory.property("FastenerPosFile").set("")
        self._factory.property("XJanchorPosFile").set("")
        self._factory.property("SJCLSanchorPosFile").set("")
        self._factory.property("SJReceiverPosFile").set("")
        self._factory.property("SJFixturePosFile").set("")

    def set_tt_name(self, tt_name):
        self._factory.property("TTName").set(tt_name)

    def mask(self):
        if not self._mask:
            self._mask = self.tool("PMTMaskConstruction")
        return self._mask

    def cd(self):
        return self.tool("DetSim1Construction")

    def wp(self):
        return self.tool("WaterPoolConstruction")

    def enable_wp_latticed(self):
        detsimalg = self.detsimalg()
        wp = self.wp()
        wp.property("enabledLatticedShell").set(True)
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

    def enable_chimney(self):
        cd = self.cd()
        cd.property("UseChimney").set(True)
    def disable_chimney(self):
        cd = self.cd()
        cd.property("UseChimney").set(False)
        self._factory.property("LowerChinmeyEnable").set(False)
        self._factory.property("TopChinmeyEnable").set(False)

    def upperchim(self):
        return self._detsimalg.createTool("UpperChimney")
    def upperchim_pos(self):
        return self._detsimalg.createTool("UpperChimneyPlacement")
    def lowerchim(self):
        return self._detsimalg.createTool("LowerChimney")

#    chim = tool(self,"Chimney")
    def enable_shutter(self):
        lowerchim=self.lowerchim()
        lowerchim.property("UseShutter").set(True)

    def set_lower_chimney(self,offsetZ,Reflec):
        lowerchim=self.lowerchim()
        lowerchim.property("UseLowerChimney").set(True)
        lowerchim.property("BlockerZ").set(offsetZ)
        lowerchim.property("Reflectivity").set(Reflec)
        cd = self.cd()
        cd.property("UseChimney").set(True)
        self._factory.property("LowerChinmeyEnable").set(True)
        self._factory.property("LowerChinmeyName").set("LowerChimney")

    def set_top_chimney(self,chim_height,Reflec):
        upperchim=self.upperchim()
        upperchim.property("UseUpperChimney").set(True)
        upperchim.property("UpperChimneyTop").set(chim_height)
        upperchim.property("Reflectivity").set(Reflec)
        upperchim_pos=self.upperchim_pos()
        upperchim_pos.property("UpperChimneyTop").set(chim_height)
        self._factory.property("TopChinmeyEnable").set(True)
        self._factory.property("TopChinmeyName").set("TopChimney")

#        print upperchim
#        print upperchim_pos

    def set_veto_pmt_mode(self, mode="CalMode"):
        self._factory.property("VetoPMTPosMode").set(mode)

    def set_qe_scale(self, scale):
        self._factory.property("ElecSimPMTQEScale").set(scale) # or 0.9

    def disable_struts_in_cd(self, t):
        # t is ["all", "strut", "fastener"]
        if t in ["all", "strut"]:
            self._factory.property("StrutPosFile").set("")
        if t in ["all", "fastener"]:
            self._factory.property("FastenerPosFile").set("")
            self._factory.property("XJanchorPosFile").set("")

    def disable_3inch_PMT(self):
        self._factory.property("3inchPMTPosFile").set("")

    def geom(self):
        if self._geom is None:
            self._geom = self.tool("GeoAnaMgr")
            if self._factory:
                self._factory.property("AnaMgrList").append("GeoAnaMgr")
        return self._geom

    def set_3inch_pmt_name(self, name):
        self._factory.property("3inchPMTName").set(name)

    def set_3inch_pmt_offset(self, value):
        self._factory.property("3inchPMTPosOffset").set(value)

    def pmtsd_mgr(self):
        if self._pmtsd_mgr is None:
            self._pmtsd_mgr = self.tool("PMTSDMgr")
        return self._pmtsd_mgr

    def optical_process(self):
        if self._optical_process is None:
            self._optical_process = self.tool("DsPhysConsOptical")
        return self._optical_process

    def em_process(self):
        if self._em_process is None:
            self._em_process = self.tool("DsG4EmLivermorePhysics")
            #self._em_process = self.tool("DsPhysConsEM")
        return self._em_process


    def enable_PMTSD_v2(self):
        # PMT SD Manager
        pmtsdmgr = self.pmtsd_mgr()
        pmtsdmgr.property("PMTSD").set("junoSD_PMT_v2")
        # Optical Physics 
        #op = self.optical_process()
        #op.property("UsePMTOpticalModel").set(False)
