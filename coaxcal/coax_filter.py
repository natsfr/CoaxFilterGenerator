#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 02:14:54 2022

@author: nats
"""
import os
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
        
    def runSim(self):
        matSpec = {'cavRad':self.cavRad * 1E3, 'wgRad': 1.45, 'coaxRad': 0.45}
        fSim = FilterSim(100e6, 700e6, 5e6, matSpec, self.coaxLPF)
        fSim.buildGeo()
        
class FilterSim:
    def __init__(self, fmin, fmax, step, spec, lpfList):
        self.unit = 1E-3
        self.FDTD = openEMS(EndCriteria=1E-4)
        self.f0 = (fmin+fmax)/2
        self.fc = (fmin+fmax)/2
        self.FDTD.SetGaussExcite(self.f0, self.fc)
        self.FDTD.SetBoundaryCond( ['MUR', 'MUR', 'MUR', 'MUR', 'MUR', 'MUR'] )
        self.CSX = CSXCAD.ContinuousStructure()
        self.FDTD.SetCSX(self.CSX)
        self.mesh = self.CSX.GetGrid()
        self.mesh.SetDeltaUnit(1e-3)
        self.mesh_res = C0/(self.f0+self.fc)/1e-3/20
        self.spec = spec
        self.lpfList = lpfList
        
    def writeXML(self):
        self.CSX.Write2XML('./geom_debug.xml')
        os.system(r'AppCSXCAD "{}"'.format('./geom_debug.xml'))
        
    def createTube(self):
        cavityRadius = self.spec['cavRad']
        extRadius = cavityRadius + 1
        z0WGRad = self.spec['wgRad']
        coaxRad = self.spec['coaxRad']
        
        # The inner radius of this shell is: rad – shell_width/2
        # The outer radius of this shell is: rad + shell_width/2
        
        gnd = self.CSX.AddMetal("GND")
        # Create end cap
        start = [0, 0, 0]
        stop = [2, 0, 0]
        
        capShellWidth = (extRadius - z0WGRad)
        capRadius = z0WGRad + capShellWidth / 2
        
        gnd.AddCylindricalShell(start, stop, capRadius, capShellWidth, priority = 2)
        start = stop
        
        lenTube = sum([i[1] for i in self.lpfList]) * 1E3
        stop = [2+lenTube, 0, 0]
        
        tubeShellWidth = (extRadius - cavityRadius)
        tubeRadius = cavityRadius + tubeShellWidth / 2
        
        gnd.AddCylindricalShell(start, stop, tubeRadius, tubeShellWidth, priority = 2)
        start = stop
        
        stop = [sum(x) for x in zip(start, [2,0,0])]
        gnd.AddCylindricalShell(start, stop, capRadius, capShellWidth, priority = 2)
        
        start = [2, 0, 0]
        stop = [2+lenTube, 0, 0]
        air = self.CSX.AddMaterial("Air", epsilon=1)
        air.AddCylinder(start, stop, cavityRadius, priority = 0)
        
        start = [0, 0, 0]
        stop = [2, 0, 0]
        ptfe = self.CSX.AddMaterial("PTFE", epsilon=2.2)
        ptfe.AddCylinder(start, stop, z0WGRad, priority = 1)
        
        start = [lenTube, 0, 0]
        stop = [lenTube+2, 0, 0]
        ptfe.AddCylinder(start, stop, z0WGRad, priority = 1)
        
        start = [0, 0, 0]
        stop = [2, 0, 0]
        coaxCopper = self.CSX.AddMetal("Coax")
        coaxCopper.AddCylinder(start, stop, coaxRad, priority = 2)
        
        start = [lenTube, 0, 0]
        stop = [lenTube+2, 0, 0]
        coaxCopper.AddCylinder(start, stop, coaxRad, priority = 2)
        return
        
    def buildGeo(self):
        self.createTube()
        self.writeXML()
        return
        