#!/usr/bin/python
#
# References:
# - https://en.wikipedia.org/wiki/B-spline#Relationship_to_piecewise/composite_B%C3%A9zier
# - https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve
#
# A series of SVG curveto (C) path commands specifies a composite bezier
# spline of cubic or 3rd order bezier curves. The endpoints of each curveto
# can have C0, C1, or C2 continuity, represented in inkscape by handles of
# arbitrary length and angle (C0), handles with a 180 deg restriction (C1), or
# handles with 180 deg and same length restriction (C2).
#
# An addintional path attribute sodipodi:nodetypes="ccszzcsc" is used to
# specify C0=c, C1=s or C2=z restrictions. There is one letter per node.
# Nodes of polygon segments (L, H, V) have letter c.
#
# Freecad sketches support higher order B-Splines. We can construct an exact
# representation of SVG curveto paths with C0, C1, or C2 properties by
# composing cubic bezier splines. The handles as known from inkscape are
# represented by lines in construction mode connecting curve points
# (two end points per cubic bezier spline) with their control points
# (two inner points per cubic bezier spline).
#
# An inkscape elliptic arc path command A has 4 groups of x,y 0 0 1 x,y
# coordinates, a total of 28.
# TODO: understand why 4 groups, conversion to one of the freecad arc types.
#       Probably easier done by evaluating the sodipodi:type,cx,xy,rx,ry,start,end,open
#
# In any case, transform attributes must be calculated, as freecad sketches
# have no transforms for their elements. Can inksvg.py resolve all that into
# proper splines?
#

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


  Example spline syntax seen in the Python-Console of FreeCaD 0.18 when drawing a cubic spline manually:

  addGeometry(Part.Circle(App.Vector(-85,192,0),App.Vector(0,0,1),10),True)    # 3
  addGeometry(Part.Circle(App.Vector(-107,160,0),App.Vector(0,0,1),10),True)   # 4
  addConstraint(Sketcher.Constraint('Radius',3,7.000000))
  addConstraint(Sketcher.Constraint('Equal',3,4))
  addGeometry(Part.Circle(App.Vector(-20,161,0),App.Vector(0,0,1),10),True)    # 5
  addConstraint(Sketcher.Constraint('Equal',3,5))
  addGeometry(Part.Circle(App.Vector(-42,193,0),App.Vector(0,0,1),10),True)    # 6
  addConstraint(Sketcher.Constraint('Equal',3,6))
  addGeometry(Part.BSplineCurve([App.Vector(-85,192),App.Vector(-107,160),App.Vector(-20,161),App.Vector(-42,193)],
              None,None,False,3,None,False),False)

  Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',3,3,7,0)
  Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',4,3,7,1)
  Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',5,3,7,2)
  Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',6,3,7,3)
  exposeInternalGeometry(7) # 7
  """

  def __init__(self, sketch, same_point, expose_int=True, circ_r=10, debug=False):
    self.ske = sketch
    self.expose_int = expose_int
    self.same_point = same_point

    self.first_point = None
    self.first_point_idx = None
    self.first_point_circ_idx = None

    self.first_circ_idx = None          # maintained by _find_or_add_circ()
    self.last_point_circ_idx = None
    self.last_point_idx = None
    self.last_point = None
    self.circ_r = circ_r
    self.debug = debug


  def line(self, p1, p2, closing=False):
    if self.same_point(p1, p2):
      return
    idx = int(self.ske.GeometryCount)
    if self.first_point is None:
      self.first_point = p1
      self.first_point_idx = (idx,1)
      self.first_point_circ_idx = None

    self.ske.addGeometry([Part.LineSegment(p1, p2)])
    if self.debug: print 'line: idx=', idx

    if self.last_point is not None and self.same_point(self.last_point, p1):
      if self.debug: print 'line: Coincident', self.last_point_idx[0], self.last_point_idx[1], idx, 1
      self.ske.addConstraint(Sketcher.Constraint('Coincident', self.last_point_idx[0], self.last_point_idx[1], idx, 1))
    if closing and self.same_point(self.first_point, p2):
      if self.debug: print 'line: Coincident Z', idx, 2, self.first_point_idx[0], self.first_point_idx[1]
      self.ske.addConstraint(Sketcher.Constraint('Coincident', idx, 2, self.first_point_idx[0], self.first_point_idx[1]))

    self.last_point = p2
    self.last_point_idx = (idx,2)
    self.last_point_circ_idx = None


  def _find_or_add_circ(self, pt, ptlist, constr=True):
    """
    When adding a constuction circle to a control point in a spline, this method
    checks the ptlist for already exising circles.
    ptlist is a list of tuples, consisting of coordinates and index.

    The index if the circle is returned.
    Caller is responsible to add newly created circles to the ptlist of subsequent invocations.
    """

    if self.debug: print "_find_or_add_circ: pt=%s, ptlist=%s" % (pt,ptlist)
    for old in ptlist:
      if self.debug: print " test %s against %s idx=%s" % (pt, old[0], old[1])
      if old is not None and old[0] is not None and old[1] is not None and self.same_point(old[0], pt):
        if self.debug: print " -> return idx=%s" % (old[1])
        return old[1]
    idx = int(self.ske.GeometryCount)
    self.ske.addGeometry(Part.Circle(pt, Normal=Base.Vector(0,0,1), Radius=self.circ_r), constr)
    if self.first_circ_idx is None:
      self.ske.addConstraint(Sketcher.Constraint('Radius', idx, self.circ_r))
      self.first_circ_idx = idx
    else:
      self.ske.addConstraint(Sketcher.Constraint('Equal',  idx, self.first_circ_idx))
    return idx


  def cubic(self, p1,h1,h2,p2, closing=False):
    """
    If closing==True, we check for hitting self.first_point
    If we are sure that no intermediate points coincide with the first_point, we can always pass closing=True.
    """

    if ((self.same_point(h1, p1) or self.same_point(h1, p2)) and
        (self.same_point(h2, p2) or self.same_point(h2, p1))):
      self.line(p1, p2, closing)
      return
    idx = int(self.ske.GeometryCount)
    self.ske.addGeometry(Part.BSplineCurve([p1, h1, h2, p2],None,None,False,3,None,False),False)
    if self.debug: print "cubic(self, p1=%s,h1=%s,h2=%s,p2=%s, closing=%s) -> idx=%s" % (p1,h1,h2,p2, closing, idx)

    # 4 circles in construction mode
    p1_circ_idx = self._find_or_add_circ(p1, [(self.last_point,self.last_point_circ_idx)])
    if self.first_point is None:
      self.first_point = p1
      self.first_point_idx = (idx,1)
      self.first_point_circ_idx = p1_circ_idx

    ptlist = [(p1,p1_circ_idx)]
    h1_circ_idx = self._find_or_add_circ(h1, ptlist)

    ptlist.append((h1,h1_circ_idx))
    h2_circ_idx = self._find_or_add_circ(h2, ptlist)

    ptlist.append((h2,h2_circ_idx))
    if closing: ptlist.append((self.first_point,self.first_point_circ_idx))
    p2_circ_idx = self._find_or_add_circ(p2, ptlist)

    conList = []
    # register the 4 circles as the 4 control points. Helps tuning...
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',p1_circ_idx,3,idx,0))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',h1_circ_idx,3,idx,1))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',h2_circ_idx,3,idx,2))
    conList.append(Sketcher.Constraint('InternalAlignment:Sketcher::BSplineControlPoint',p2_circ_idx,3,idx,3))
    self.ske.addConstraint(conList)

    if self.expose_int:
      self.ske.exposeInternalGeometry(idx)

    if self.debug: print "cubic: test self.same_point(self.last_point=%s, p1=%s)" % (self.last_point, p1)
    if self.last_point is not None and self.same_point(self.last_point, p1):
      if self.debug: print "cubic: Coincident", self.last_point_idx[0], self.last_point_idx[1], idx, 1
      if self.last_point_circ_idx == p1_circ_idx:
        # Be very careful with duplicate constraints. FreeCAD does not like that at all!
        if self.debug: print " -> skipped, they coincide already, as they share a circle"
        pass
      else:
        self.ske.addConstraint(Sketcher.Constraint('Coincident', self.last_point_idx[0], self.last_point_idx[1], idx, 1))
    if closing and self.same_point(self.first_point, p2):
      if self.debug: print 'cubic: Coincident Z', idx, 2, self.first_point_idx[0], self.first_point_idx[1]
      if self.first_point_circ_idx == p2_circ_idx:
        # Be very careful with duplicate constraints. FreeCAD does not like that at all!
        if self.debug: print " -> skipped, they coincide already, as they share a circle"
      else:
        self.ske.addConstraint(Sketcher.Constraint('Coincident', idx, 2, self.first_point_idx[0], self.first_point_idx[1]))

    self.last_point = p2
    self.last_point_idx = (idx,2)
    self.last_point_circ_idx = p2_circ_idx


def vec(p): return Base.Vector(p[0], -p[1], 0)

def _same_point(p1, p2, eps=0.001):
  if abs(p1[0]-p2[0]) > eps: return False
  if abs(p1[1]-p2[1]) > eps: return False
  return True

spt = SubPathTracker(ske, _same_point, debug=True)
# first part "m 0,10 c 60,0 -20,50  40,50"
p1, h1, h2, p2 = ( (0,10), (60,10), (-20,50), (40,50) )
spt.cubic(vec(p1),vec(h1),vec(h2),vec(p2))

# second part "m 40,50 c 40,0 0,-40 40,-40"
p1, h1, h2, p2 = ( (40,50), (80,50), (40,10), (80,10) )
spt.cubic(vec(p1),vec(h1),vec(h2),vec(p2), closing=True)

p1, p2 = ( (80,10), (40,-30) )
spt.line(vec(p1),vec(p2), closing=True)

p1, h1, h2, p2 = ( (40,-30), (20,-50), (-20,10), (0,10) )
spt.cubic(vec(p1),vec(h1),vec(h2),vec(p2), closing=True)

doc.saveAs('./sketch_spline_py.fcstd')
print('done ./sketch_spline_py.fcstd')


