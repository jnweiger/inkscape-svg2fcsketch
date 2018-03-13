#!/usr/bin/python
 
import os
import sys
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')
 
from FreeCAD import Base
import Part
 
# # Parametric Modelling
doc = FreeCAD.newDocument('box')
 
body  = Part.makeBox(2,2,2)             # Create a Box
transvec = Base.Vector(20,0,0)          # Translation Vector
body.translate(transvec)                # Translate the body
Part.show(body)                         # Add the body to ActiveDocument
 
doc.saveAs('./box-py.fcstd')

