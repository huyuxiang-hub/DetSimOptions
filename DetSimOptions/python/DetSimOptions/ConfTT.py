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

        self._gdml = None

    def configure(self):
        # == G4Svc ==
        self._g4svc = self._task.createSvc("G4Svc")

        # == DetSimOptions ==
        self._factory = self._task.createSvc("DetSim0Svc")
        self._factory.property("AnaMgrList").set(["DataModelWriter",
                                                  "GenEvtInfoAnaMgr",
                                                  "NormalAnaMgr", 
                                                  "DepositEnergyAnaMgr",
                                                  "DepositEnergyTTAnaMgr"])
        self._factory.property("CDName").set("DetSim1")
        #self._factory.property("PMTPosFile").set(DetSimOptions.data_load("PMTPos_Acrylic.csv"))
        #self._factory.property("StrutPosFile").set(DetSimOptions.data_load("Strut_Acrylic.csv"))
        #self._factory.property("FastenerPosFile").set(DetSimOptions.data_load("Strut_Acrylic.csv"))

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

        return self._detsimalg.createTool(toolname)

    def set_gdml_output(self, filename):
        # If already setup, just return
        if self._gdml and self._gdml != filename:
            raise Exception("You have already set the gdml output.")
        self._gdml = filename
        # if already exists, skip this
        if os.path.exists(self._gdml):
            print "The %s already exists, skip generate" %self._gdml
            return
        # configure the GeoAnaMgr 
        geo = self.tool("GeoAnaMgr")
        geo.property("GdmlEnable").set(True)
        geo.property("GdmlOutput").set(self._gdml)
        if self._factory:
            self._factory.property("AnaMgrList").append("GeoAnaMgr")

    def mask(self):
        if not self._mask:
            self._mask = self.tool("PMTMaskConstruction")
        return self._mask

    def set_tt_name(self, ttname):
        if not self._detsimalg:
            self.configure()
        self._factory.property("TTName").set(ttname)
