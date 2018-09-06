#!/usr/bin/python
#
#  https://en.wikipedia.org/wiki/B-spline#Relationship_to_piecewise/composite_B%C3%A9zier
 
import os
import sys
sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')
 
from FreeCAD import Base
import Part, Sketcher   # causes SEGV if Base is not yet imported from FreeCAD
# import ProfileLib.RegularPolygon as Poly
doc = FreeCAD.newDocument('sketch_spline')
ske = doc.addObject('Sketcher::SketchObject', 'Sketch22')

# d="m 0,-10 c 60,0 -20,40 40,40 40,0 0,-40 40,-40"
# inkscape:connector-curvature="0"
# sodipodi:nodetypes="ccc"

class SubPathTracker():
  """
  Track status of a subpath. In SVG, a path consists of disconnected subpaths.
  Each subpath consists of a connected list of straight lines or cubic spline segments.

  A SubPathTracker object remembers the first and last last point of a subpath.
  When adding a line() or cubic(), it continues the current subpath 
  assuming same_point(last_point, p1) == True, otherwise it is an error.

  Control points of splines have circles in construction mode on them. All the same size.
  End points of lines don't have these, unless they are also control points of adjacent splines.
  
  If any two or three of the control points coincide, duplicate circles are avoided, an coincidence restrictions are added instead.
  We never convert a cubic spline to a quadratic spline, as they are slightly different in shape.
  
  Adjacent ends of segments have coincidence restrictions on them, no matter if the segments are of type line or type cubic.
  Duplicate circles on adjacent ends of two cubics are avoided.

  If the very last point of a subpath coincides with the very first point, a coincidence restriction is added, and duplicate circles 
  are again avoided.

  Other coincidences within the subpath are ignored and may prodice duplicate circles (serving their own purpose each).
  
  """
  def __init__(self, sketch, same_point, expose_int=True, circ_r=10):
    self.ske = sketch
    self.expose_int = expose_int
    self.same_point = same_point
    self.first_circ_idx = None
    self.last_circ_idx = None
    self.circ_r = circ_r

  def _round_sigdigs(val, n=2):
    """
    Rounds the value to at most n significant digits.
    0.00426221  -> 0.0043
    3.78        -> 3.8
    997         -> 1000
    994         -> 990
    """
    n = int(max(1, n))
    exp10 = math.floor(math.log10(abs(val))) if val != 0 else 0
    decimal_shifter = float(10**exp10)
    val = round(val/decimal_shifter, n-1)
    # this final round is mathematically useless, but
    # avoids _round_sigdigs(0.000491) -> 0.0004900000000000001
    return round(val*decimal_shifter, -exp10+n-1)       

  def _average_handle_length(self, sp):
    tot = 0
    cnt = 0
    for tri in sp:
      (h1, p, h2) = tri
      if not self.same_point(h1, p):
        tot += math.sqrt( (h1[0]-p[0])*(h1[0]-p[0]) + (h1[1]-p[1])*(h1[1]-p[1]) )
        cnt += 1
      if not self.same_point(h2, p):
        tot += math.sqrt( (h2[0]-p[0])*(h2[0]-p[0]) + (h2[1]-p[1])*(h2[1]-p[1]) )
        cnt += 1
    if (cnt > 0 and tot > 0):
      return self._round_sigdigs(tot/cnt, 2)
    return 10


  def line(self, p1, p2):
    return if self.same_point(p1, p2)
    idx = int(self.ske.GeometryCount)
    self.ske.addGeometry([Part.LineSegment(p1, p2)])
    if self.first_point is None:
      self.first_point = p1
      self.first_point_idx = (idx,0)
    self.last_circ = None
    self.last_point = p2
    self.last_point_idx = (idx,1)


  def cubic(self, p1,h1,h2,p2):
    if self._same_point(h1, p1) and self._same_point(h2, p2):
      self.line(p1, p2):
      return
    idx = int(self.ske.GeometryCount)
    # 4 circles in construction mode with equal diameter
    self.ske.addGeometry(Part.Circle(p1, Normal=Base.Vector(0,0,1), Radius=self.circ_r),True)
    self.ske.addGeometry(Part.Circle(h1, Normal=Base.Vector(0,0,1), Radius=self.circ_r),True)
    self.ske.addGeometry(Part.Circle(h2, Normal=Base.Vector(0,0,1), Radius=self.circ_r),True)
    self.ske.addGeometry(Part.Circle(p2, Normal=Base.Vector(0,0,1), Radius=self.circ_r),True)
    self.ske.addConstraint(Sketcher.Constraint('Radius', idx+0, self.circ_r))
    self.ske.addConstraint(Sketcher.Constraint('Equal',  idx+0, idx+1))
    self.ske.addConstraint(Sketcher.Constraint('Equal',  idx+0, idx+2))
    self.ske.addConstraint(Sketcher.Constraint('Equal',  idx+0, idx+3))
  
    self.ske.addGeometry(Part.BSplineCurve([p1, h1, h2, p2],None,None,False,3,None,False),False)
    conList = []
    # register the 4 circles as the 4 control points. Helps tuning...
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',idx+0,3,idx+4,0))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',idx+1,3,idx+4,1))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',idx+2,3,idx+4,2))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',idx+3,3,idx+4,3))
    self.ske.addConstraint(conList)
    if self.expose_int:
      self.ske.exposeInternalGeometry(idx+4)


def vec(p): return Base.Vector(p[0], -p[1], 0)
sp = SplineTracker(ske)
# first part "m 0,10 c 60,0 -20,50  40,50"
p1, h1, h2, p2 = ( (0,10), (60,10), (-20,50), (40,50) )
sp.cubic(vec(p1),vec(h1),vec(h2),vec(p2))

# second part "m 40,50 c 40,0 0,-40 40,-40"
p1, h1, h2, p2 = ( (40,50), (80,50), (40,10), (80,10) )
sp.cubic(vec(p1),vec(h1),vec(h2),vec(p2))
 
doc.saveAs('./sketch_spline_py.fcstd')
print('done ./sketch_spline_py.fcstd')


