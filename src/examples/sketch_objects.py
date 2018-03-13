#!/usr/bin/python
# FROM: 
#  https://github.com/FreeCAD/FreeCAD/blob/master/src/Mod/Sketcher/TestSketcherApp.py
#  /usr/lib/freecad-daily/Mod/Sketcher/SketcherExample.py
 
import os, sys, math, re
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')

outdir = '.'
docname = 'sketch_objects'
if len(sys.argv) > 1: 
  name = sys.argv[1]
docname = re.sub('\.fcstd$', '', docname)
m = re.match('^(.*)/(.*?)$', docname)
if m: 
  outdir = m.group(1)
  docname = m.group(2)

 
from FreeCAD import Base
import Part, Sketcher
import ProfileLib.RegularPolygon as Poly
 
doc = FreeCAD.newDocument(docname)
 
ske = doc.addObject('Sketcher::SketchObject', 'Sketch')         
# FIXME: sketch is rendered invisible
# We need a DocumentViewProviderObject to set Visibility=True
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
    #Part.ArcOfCircle(Part.Circle(Base.Vector(4,7,0), Base.Vector(0,0,1), 1), math.pi/2, -math.pi)
    #
    #
    # Long axis SW: (N, -N, 0)
    # Short axis NE: (N, N, 0)
    # StartPoint North: 3/4 pi  (counted CCW from long axis)
    # EndPoint   NNE:   -9/8 pi (counted CW from long axis)
    # FIXME: the angles are distroted according to the tangents of the contact point.
    # North is more like 5/8 pi but a bit less.
    Part.ArcOfEllipse(Part.Ellipse(App.Vector(8,-8,0),App.Vector(3,3,0),App.Vector(0,0,0)),-math.pi*(8+1.)/8, math.pi*(5./8))
]
ske.addGeometry(geo, False)     # False: normal, True: blue construction
ske.exposeInternalGeometry(i+5) # long and short axis, focal points of Ellipse

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

#ske.addConstraint(con)
doc.recompute()

file = outdir+'/'+docname+'.fcstd'
doc.saveAs(outdir+'/'+docname+'.fcstd')
print("%s written." % file)

