#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: lintao

import os.path
import Sniper
import DetSimOptions
from .ConfAcrylic import ConfAcrylic

class ConfBalloon(ConfAcrylic):

    def __init__(self, task):
        super(ConfBalloon, self).__init__(task)

    def configure(self):
        # == G4Svc ==
        self._g4svc = self._task.createSvc("G4Svc")

        # == DetSimOptions ==
        self._factory = self._task.createSvc("DetSim0Svc")
        self._factory.property("AnaMgrList").set(["DataModelWriter",
                                                  "GenEvtInfoAnaMgr", 
                                                  "NormalAnaMgr",
                                                  "DepositEnergyAnaMgr"])
        self._factory.property("CDName").set("DetSim2")
        self._factory.property("PMTPosFile").set(DetSimOptions.data_load("PMTPos_Balloon.csv"))
        self._factory.property("3inchPMTPosFile").set(DetSimOptions.data_load("3inch_pos.csv"))
        self._factory.property("StrutPosFile").set(DetSimOptions.data_load("Strut_Balloon.csv"))

        # == DetSimAlg ==
        self._detsimalg = self._task.createAlg("DetSimAlg")
        self._detsimalg.property("DetFactory").set(self._factory.objName())
        self._detsimalg.property("RunMac").set("run.mac")
        
