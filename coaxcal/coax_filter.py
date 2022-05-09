#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 02:14:54 2022

@author: nats
"""

import numpy as np

from CSXCAD import CSXCAD

from openEMS.openEMS import openEMS
from openEMS.physical_constants import *

class CoaxFilter:
    def __init__(self, maxRad, minRad, cavRad, Er, Ur):
        self.maxRad = maxRad
        self.minRad = minRad
        self.cavRad = cavRad
        self.Er = Er
        self.E0 = 8.854187E-12
        self.Ur = Ur
        self.U0 = 4E-7*np.pi
        self.lpfList = []
        self.coaxLPF = []
        self.capUnit = 1E12
        self.capStrUnit = "pF"
        self.indUnit = 1E9
        self.indStrUnit = "nH"
        
    def getCap(self, length):
        return 2*np.pi*self.Er*self.E0/(np.log(self.cavRad/self.maxRad)) * length
    
    def getInd(self, length):
        return (self.Ur*self.U0)/(2*np.pi)*np.log(self.cavRad/self.minRad) * length
    
    def getCapLen(self, cap):
        return (cap * np.log(self.cavRad/self.maxRad))/(2*np.pi*self.Er*self.E0)

    def getIndLen(self, ind):
        return (ind * 2 * np.pi)/(self.Ur*self.U0*np.log(self.cavRad/self.minRad))
    
    def setLPF(self, lpfList):
        self.lpfList = lpfList
        self.order = len(lpfList)
        
    def strictConvert(self):
        return
    
    def setCoaxLPF(self, coaxLPF):
        self.coaxLPF = coaxLPF
        
class FilterSim:
    def __init__(self, fmin, fmax, step, spec, lpfList):
        self.unit = 1E-3
        self.FDTD = openEMS(NrTs=30e3, EndCriteria=1E-4)
        self.f0 = (fmin+fmax)/2
        self.fc = (fmin+fmax)/2
        self.FDTD.SetGaussExcite(f0, fc)
        self.FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )
        self.CSX = CSXCAD.ContinuousStructure()
        self.FDTD.SetCSX(self.CSX)
        self.mesh = CSX.GetGrid()
        self.mesh.SetDeltaUnit(1e-3)
        self.mesh_res = C0/(f0+fc)/1e-3/20
        self.spec = spec
        self.lpfList = lpfList
        
    def createTube(self):
        cavityRadius = self.spec['cavRad']
        extRadius = cavityRadius + 1
        z0WGRad = self.spec['wgRad']
        
        gnd = self.CSX.AddMetal("GND")
        # Create end cap
        start = [0, 0, 0]
        stop = [2, 0, 0]
        gnd.AddCylindricalShell(10, start, stop, extRadius, extRadius-z0WGRad)
        start = stop
        lenTube = sum([i[1] for i in self.lpfList])
        stop = [2+lenTube, 0, 0]
        gnd.AddCylindricalShell(10, start, stop, extRadius, extRadius-z0WGRad)
        return
        
    def buildGeo(self, lpfList):
        return
        