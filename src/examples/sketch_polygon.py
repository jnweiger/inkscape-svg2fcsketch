#!/usr/bin/python
 
import os
import sys
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')
 
from FreeCAD import Base
import ProfileLib.RegularPolygon as Poly
 
doc = FreeCAD.newDocument('sketch_poly')
 
ske = doc.addObject('Sketcher::SketchObject', 'Sketch')
# ske.Placement = Base.Placement(Base.Vector(0, 0, 0),
#                               Base.Rotation(-0.707107, 0, 0, -0.707107))
Poly.makeRegularPolygon('Sketch', sides=5,
                        centerPoint=Base.Vector(0, 0, 0),
                        firstCornerPoint=Base.Vector(10,0,0),
                        construction=False)

doc.saveAs('./sketch_poly-py.fcstd')

