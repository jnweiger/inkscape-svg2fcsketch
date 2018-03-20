#!/usr/bin/python
 
import os, math
import sys
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')

from FreeCAD import Base
import Part, Sketcher
 
doc = FreeCAD.newDocument('sketch_rect')
 
ske = doc.addObject('Sketcher::SketchObject', 'Sketch')
# ske.Placement = Base.Placement(Base.Vector(0, 0, 0),
#                                Base.Rotation(-0.707107, 0, 0, -0.707107))

# Base.Rotation(Base.Vector(0,0,1),45).multVec(Base.Vector(0,10,0))
# ->  Vector (-7.0710678118654755, 7.071067811865475, 0.0)

## Rotation 15deg:
# a = 15/180.*math.pi
# m = Base.Matrix(math.cos(-a), -math.sin(-a), 0, 0, math.sin(-a), math.cos(-a), 0, 0)
# v = m.multiply(Base.Vector(10, 0, 0))
# v = Vector (9.659258262890683, -2.5881904510252074, 0.0)


def matrix2d(svgmat):
  """
  Convert a 2D matrix from SVG into a FreeCAD 3D Matrix.
  e.g. mat = [[0.9659258262890683, 0.25881904510252074, 0.0], 
              [-0.25881904510252074, 0.9659258262890683, 0.0]]
  results in Matrix ((0.965926,0.258819,0,0),(-0.258819,0.965926,0,0),(0,0,1,0),(0,0,0,1))
  """
  # FIXME: is the handling of svgmat[*][2] correct?
  return Base.Matrix(svgmat[0][0], svgmat[0][1], 0, svgmat[0][2], 
                     svgmat[1][0], svgmat[1][1], 0, svgmat[1][2])

def vec2d(x, y, m):
  """
  Converts a 2D Vector from SVG into a FreeCAD 3D vector applying Base.Matrix m.
  """
  return m.multiply(Base.Vector(x, y, 0))


(x, y, w, h) = (39.276257, 86.713783, 61.232143, 43.845234)
# svg style matrix, with a 15deg rotation:
m = matrix2d([[0.9659258262890683, 0.25881904510252074, 0.0], [-0.25881904510252074, 0.9659258262890683, 0.0]])
m.move(0,50,0)  # inplace: move origin from top of SVG page to bottom of SVG page.
m.scale(1,-1,1)  # inplace: flip the Y-Axis, SVG has left-handed coordinate system.

i = int(ske.GeometryCount)
ske.addGeometry([
    Part.LineSegment(vec2d(x  ,y  ,m), vec2d(x+w,y  ,m)),
    Part.LineSegment(vec2d(x+w,y  ,m), vec2d(x+w,y+h,m)),
    Part.LineSegment(vec2d(x+w,y+h,m), vec2d(x  ,y+h,m)),
    Part.LineSegment(vec2d(x  ,y+h,m), vec2d(x  ,y  ,m))
], False)
ske.addConstraint([
    Sketcher.Constraint('Coincident',i+0,2, i+1,1),
    Sketcher.Constraint('Coincident',i+1,2, i+2,1),
    Sketcher.Constraint('Coincident',i+2,2, i+3,1),
    Sketcher.Constraint('Coincident',i+3,2, i+0,1),
])

doc.saveAs('./sketch_rect.fcstd')

