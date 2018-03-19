#!/usr/bin/python
# FROM: 
#  https://github.com/FreeCAD/FreeCAD/blob/master/src/Mod/Sketcher/TestSketcherApp.py
#  /usr/lib/freecad-daily/Mod/Sketcher/SketcherExample.py
 
import os, sys, math
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')
 
from FreeCAD import Base
import Part, Sketcher
import ProfileLib.RegularPolygon as Poly
 
doc = FreeCAD.newDocument('sketch_lines')
 
ske = doc.addObject('Sketcher::SketchObject', 'Sketch')
ske.Placement = Base.Placement(Base.Vector(0, 0, 0),
                               Base.Rotation(0, 0, 0, 1))
i = int(ske.GeometryCount)      # 0
geo = [
    Part.LineSegment(Base.Vector(4,8,0),Base.Vector(9,8,0)),
    Part.LineSegment(Base.Vector(9,8,0),Base.Vector(9,2,0)),
    Part.LineSegment(Base.Vector(9,2,0),Base.Vector(3,2,0)),
    Part.LineSegment(Base.Vector(3,2,0),Base.Vector(3,7,0)),
    Part.Circle(Center=Base.Vector(6,5,0), Normal=Base.Vector(0,0,1), Radius=2),
    # North: math.pi/2
    # SouthEast: -math.pi/4
    # West: math.pi, -math.pi
    # East: 0
    Part.ArcOfCircle(Part.Circle(Base.Vector(4,7,0), Base.Vector(0,0,1), 1), math.pi/2, -math.pi)
]
ske.addGeometry(geo, False)
print("GeometryCount changed from %d to %d" % (i, int(ske.GeometryCount)))

con = [
    Sketcher.Constraint('Coincident',i+0,2, i+1,1),
    Sketcher.Constraint('Coincident',i+1,2, i+2,1),
    Sketcher.Constraint('Coincident',i+2,2, i+3,1),
    Sketcher.Constraint('Coincident',i+3,2, i+5,2),
]

#ske.addConstraint(Sketcher.Constraint('Horizontal',0))
#ske.addConstraint(Sketcher.Constraint('Horizontal',2))
#ske.addConstraint(Sketcher.Constraint('Vertical',1))
#ske.addConstraint(Sketcher.Constraint('Vertical',3))

ske.addConstraint(con)
#doc.recompute()

doc.saveAs('./sketch_lines-py.fcstd')

