#! /usr/bin/python
#
# (C) 2018 juergen@fabmail.org
# Distribute under GPL-2.0 or ask.
#
# References:
#  https://github.com/FreeCAD/FreeCAD/blob/master/src/Mod/Sketcher/TestSketcherApp.py
#  /usr/lib/freecad-daily/Mod/Sketcher/SketcherExample.py
#  https://en.wikipedia.org/wiki/Rytz%27s_construction#Computer_aided_solution
#  http://wiki.inkscape.org/wiki/index.php/Python_modules_for_extensions
#  https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve
#  https://en.wikipedia.org/wiki/B-spline#Relationship_to_piecewise/composite_B%C3%A9zier
#
# v0.1 jw, initial draft refactoring inksvg to make it fit here.
# v0.2 jw, Introducing class SketchPathGen to seperate the sketch generator from the svg parser.
# v0.3 jw, correct _coord_from_svg() size and offset handling. Suppress
#          silly version printing, that would ruin an inkscape extension.
# V0.4 jw, Added GuiDocument.xml for visibility and camera defaults.
#          Using BoundBox() to compute camera placement.
# V0.5 jw, objEllipse() done correctly with _ellipse_vertices2d()
# V0.6 jw, objArc() done. ArcOfCircle() is a strange beast with rotation and mirroring.
# V0.7 jw, pathString() done.
# V0.8 jw, imported class SubPathTracker() from src/examples/sketch_spline.py
#

from optparse import OptionParser
import os, sys, math, re

sys_platform = sys.platform.lower()
if sys_platform.startswith('win'):
  sys.path.append('C:\Program Files\Inkscape\share\extensions')
elif sys_platform.startswith('darwin'):
  sys.path.append('~/.config/inkscape/extensions')
else:   # Linux
  sys.path.append('/usr/share/inkscape/extensions/')

import cubicsuperpath

sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')

verbose=-1       # -1=quiet, 0=normal, 1=babble
epsilon = 0.00001

if verbose <= 0:
  os.dup2(1,99)         # hack to avoid silly version string printing.
  f = open("/dev/null", "w")
  os.dup2(f.fileno(), 1)

# The version printing code has
# src/App/Application.cpp: if (!(mConfig["Verbose"] == "Strict"))
# but we cannot call SetConfig('Verbose', 'Strict') early enough.
from FreeCAD import Base, BoundBox
sys.stdout.flush()      # push silly version string into /dev/null

if verbose <= 0:
  f.close()
  os.dup2(99,1)         # back in cansas.

import Part, Sketcher   # causes SEGV if Base is not yet imported from FreeCAD
import ProfileLib.RegularPolygon as Poly

# CAUTION: Keep in sync with with svg2fcsketch.inx ca. line 3 and line 24
__version__ = '0.8'

## INLINE_BLOCK_START
# for easier distribution, our Makefile can inline these imports
# sys.path.append('/home/src/github/jnweiger/inkscape-thunderlaser/src/')
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/../inksvg/src/')
from inksvg import InkSvg, PathGenerator
## INLINE_BLOCK_END

parser = OptionParser(usage="\n    %prog [options] SVGFILE [OUTFILE]\n\nTry --help for details.")
parser.add_option("-o", "--outfile", dest="outfile",
                         help="write fcstd to OUTPUT. Default: stdout (unless it is a tty)", metavar="OUTPUT")
parser.add_option("-i", "--id", "--ids", dest="ids",
                        action="append", type="string", default=[],
                        help="Select svg object(s) by id attribute. Use multiple times or combine with comma. Default: root object, aka all")
parser.add_option("--selected-nodes", dest="selected_nodes",
                        action="append", type="string", default=[],
                        help="id:subpath:position of selected nodes, if any")
                        # TODO: check if inkscape is really passing us these
parser.add_option("-V", "--version",
                         action="store_true", dest="version", default=False,
                         help="Print version numbers only.")
parser.add_option("-v", "--verbose",
                         action="store_true", dest="verbose", default=False,
                         help="Be verbose. Default: silent.")
parser.add_option("-e", "--expose-internal-geometry",
                         action="store_true", dest="expose_internal_geometry", default=False,
                         help="Expose internal geometry for Splines and Ellipses. Default: False.")
(options, args) = parser.parse_args()

if options.version:
  print("InkSvg %s\n%s %s\n" % (InkSvg.__version__, sys.argv[0], __version__))
  sys.exit(0)
if len(args) < 1:
  parser.error("Input svg file missing.")

svgfile = args[0]

outdir = '.'
if len(args) > 1:
  fcstdfile = args[1]
else:
  if options.outfile:
    fcstdfile = options.outfile
  else:
    fcstdfile = re.sub('\.svg$', '.fcstd', svgfile, re.I)
docname = re.sub('\.fcstd$', '', fcstdfile, re.I)
docname = re.sub('^.*/', '', docname)

if not options.outfile:
  if sys.stdout.isatty():
    print "ERROR: stdout isatty. Please use option --outfile or redirect stdout."
    sys.exit(1)


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


class SketchPathGen(PathGenerator):
  """
  Generate XML code for a FreeCAD sketch.
  """
  def __init__(self, ske, yflip=True, expose_int=True):
    self.ske = ske
    self.m = None
    self.yflip = yflip
    self.expose_int = expose_int
    self.bbox = BoundBox()

    # prepare counters how many objects of each type we generate.
    self.stats = {}
    for s in ('pathList', 'pathString', 'objRect', 'objRoundedRect',
              'objArc', 'objEllipse'):
      self.stats[s] = 0


  def _coord_from_svg(self, m=None):
    """
    Use SVG properties set by InkSVG.handleViewBox()
    to define the transformation matrix from dpi'ish SVG space
    into metric FreeCAD space.
    """
    xs = 1.0
    if self._svg.docTransform is not None:
      xs = self._svg.docTransform[0][0]
      if abs(xs) < epsilon:
        xs = 1.0        # avoid divison by zero.
    ys = xs
    if self.yflip: ys = -ys

    if m is None: m = self.m
    m.move(0, -self._svg.docHeight, 0)
    m.scale(1/xs, 1/ys, 1)
    return m

  def _ellipse_vertices2d(self, C, P, Q):
    """
    Compute two vertices of the ellipse, major and minor axis.
    Given two conjugated half diameters P, Q.
    Vectors are expected as FreeCAD Base.Vector() assuming that all Z=0.

    Using a computation derived from Rytz's construction as
    seen in https://en.wikipedia.org/wiki/Rytz%27s_construction.
    Returns (V1, V2, V3, V4) where V1-V4 is the diameter of the major axis
    and V2-V3 is the diameter of the minor axis.
    """
    f0 = C
    f1 = P-C
    f2 = Q-C
    # det = (f1*f1 - f2*f2) / ( 2 * f1 * f2 )     # raise NaN, if flat.
    idet = (2*f1*f2) / (f1*f1 - f2*f2) # raise NaN, if flat or colinear.

    def math_cot(a): return 1/math.tan(a)

    def p(t): return f0 + f1 * math.cos(t) + f2 * math.sin(t)

    # cot(2*t0) = det
    # tan(2*t0) = 1/det
    # tan(2*t0) = idet
    # 2*t0 = atan(idet)
    t0 = 0.5* math.atan(idet)
    V1 = p(t0)
    V2 = p(t0 + 0.5*math.pi)
    V3 = p(t0 - 0.5*math.pi)
    V4 = p(t0 + math.pi)
    if (V1-V4).Length < (V2-V3).Length:
      # V1-V4 should be the major axis.
      (V1,V4, V2,V3) = (V2,V3, V1,V4)
    return (V1, V2, V3, V4)


  def _decompose_matrix2d(self, m=None):
    """
    Decompose a 4x4 matrix into 2d translation vector,
    2d scale vector, and rotation angle.

    Inspired by https://math.stackexchange.com/questions/13150/extracting-rotation-scale-values-from-2d-transformation-matrix
    """
    if m is None: m = self.m
    a =  m.A11
    b =  m.A12
    xc = m.A14          # FIXME: check if correct in _matrix_from_svg()
    c =  m.A21
    d =  m.A22
    yc = m.A24          # FIXME: check if correct in _matrix_from_svg()
    sign_a = -1 if a < 0 else 1
    sign_d = -1 if d < 0 else 1

    sx = sign_a * math.sqrt(a*a + b*b)
    sy = sign_d * math.sqrt(c*c + d*d)

    sign = math.atan(-c / a)
    rad  = math.acos(a / sx)
    deg  = rad * 180. / math.pi
    if (deg > 90 and sign > 0) or (deg < 90 and sign < 0):
      deg = 360 - deg
    rad = deg / 180. * math.pi

    return (Base.Vector(xc, yc), Base.Vector(sx, sy), (rad, deg))


  def _matrix_from_svg(self, svgmat, coordcvt=True):
    """
    Convert a 2D matrix from SVG into a FreeCAD 3D Matrix.
    e.g. mat = [[0.9659258262890683, 0.25881904510252074, 0.0],
               [-0.25881904510252074, 0.9659258262890683, 0.0]]

    If coordcvt is True, then coordinate system conversion from
    SVG to FreeCAD is applied to the matrix. Otherwise only datatype
    conversion is performed.

    Returns:
    e.g. Matrix ((0.965926,0.258819,0,0),(-0.258819,0.965926,0,0),(0,0,1,0),(0,0,0,1))
    """
    # FIXME: is the handling of svgmat[*][2] correct?
    self.m = Base.Matrix(svgmat[0][0], svgmat[0][1], 0, svgmat[0][2],
                         svgmat[1][0], svgmat[1][1], 0, svgmat[1][2])

    if coordcvt:
      self._coord_from_svg()
    return self.m


  def _from_svg(self, x, y, m=None, bbox=True):
    """
    Converts a 2D Vector from SVG into a FreeCAD 3D vector applying Base.Matrix m.
    """
    if m is None: m = self.m
    v = m.multiply(Base.Vector(x, y, 0))
    if bbox: self.bbox.add(v)
    return v

  def _same_point(self, p1, p2, eps=epsilon):
    if p1 is None or p2 is None: return True
    if abs(p1[0]-p2[0]) > eps: return False
    if abs(p1[1]-p2[1]) > eps: return False
    return True

  def _round_sigdigs(self, val, n=2):
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
    return round(val*decimal_shifter, int(-exp10+n-1))

  def _average_handle_length(self, sp):
    (tra,sca,rot) = self._decompose_matrix2d()
    sca = 0.5 * (abs(sca[0]) + abs(sca[1]))

    tot = 0
    cnt = 0
    for tri in sp:
      (h1, p, h2) = tri
      if not self._same_point(h1, p):
        tot += math.sqrt( (h1[0]-p[0])*(h1[0]-p[0]) + (h1[1]-p[1])*(h1[1]-p[1]) )
        cnt += 1
      if not self._same_point(h2, p):
        tot += math.sqrt( (h2[0]-p[0])*(h2[0]-p[0]) + (h2[1]-p[1])*(h2[1]-p[1]) )
        cnt += 1
    if (cnt > 0 and tot > 0):
      return self._round_sigdigs(sca*tot/cnt, 2)
    return self._round_sigdigs(sca*10, 2)


  def pathString(self, d, node, mat):
    """
    d is expected formatted as an svg path string here.
    d = "M 30.994048,129.93452 72.571427,88.357143 V 129.93452 H 127" means
    path = [
      [
        #       handle_0                 point                handle_1
        [[30.994048, 129.93452], [30.994048, 129.93452], [30.994048, 129.93452]],
        [[72.571427, 88.357143], [72.571427, 88.357143], [72.571427, 88.357143]],
        [[72.571427, 129.93452], [72.571427, 129.93452], [72.571427, 129.93452]],
        [[127.0, 129.93452],     [127.0, 129.93452],     [127.0, 129.93452]]
      ]
    ]
    """
    self._matrix_from_svg(mat)
    path = cubicsuperpath.parsePath(d)
    for subpath in path:
      spt = SubPathTracker(self.ske, lambda a,b: self._same_point(a, b), self.expose_int,
                           circ_r=0.1*self._average_handle_length(subpath), debug=(verbose > 0))

      # These are the off by one's: four points -> three lines -> two constraints.
      j = 0
      while j < len(subpath)-1:
        # lists of three, http://wiki.inkscape.org/wiki/index.php/Python_modules_for_extensions#cubicsuperpath.py
        (h0,p1,h1) = subpath[j]
        j = j+1
        while j < len(subpath):
          (h2,p2,h3) = subpath[j]
          if not self._same_point(p1, p2):
            break               # no null-segments, please!
          j += 1
        if j >= len(subpath):
          break                 # nothing left.

        spt.cubic( self._from_svg(p1[0], p1[1]), self._from_svg(h1[0], h1[1]),
                   self._from_svg(h2[0], h2[1]), self._from_svg(p2[0], p2[1]), closing=(j+1 == len(subpath)) )

      self.stats['pathString'] += 1     # count subpaths


  def pathList(self, d, node, mat):
    """
    d is expected as an [[cmd, [args]], ...] arrray
    """
    print(d, node, mat)
    i = int(self.ske.GeometryCount)      # 0
    geo = []
    geo.append(Part.LineSegment(Base.Vector(4,8,0),Base.Vector(9,8,0)))
    geo.append(Part.LineSegment(Base.Vector(9,8,0),Base.Vector(9,2,0)))
    self.ske.addGeometry(geo, False)
    print("GeometryCount changed from %d to %d" % (i, int(self.ske.GeometryCount)))
    print("not impl. simplePath: ", d, node, mat)


  def objRoundedRect(self, x, y, w, h, rx, ry, node, mat):
    """
    Construct four arcs, one for each corner, and
    connect them with line segments, if space permits.
    Connect them directly otherwise.
    """
    if rx == 0: rx = ry
    if ry == 0: ry = rx
    if rx < epsilon or ry < epsilon:
      return self.objRect(x, y, w, h, node, mat)
    if 2*rx > w-epsilon: rx = 0.5*(w-epsilon)   # avoid Part.OCCError: Both points are equal" on LineSegment #12
    if 2*ry > h-epsilon: ry = 0.5*(h-epsilon)

    if verbose > 0: print("objRoundedRect: ", x, y, w, h, rx, ry, node.get('id'), mat)

    self._matrix_from_svg(mat)
    i = int(self.ske.GeometryCount)
    ske.addGeometry([
      # construction outline of the box
      Part.LineSegment(self._from_svg(x  ,y  ), self._from_svg(x+w,y  )),               # 0
      Part.LineSegment(self._from_svg(x+w,y  ), self._from_svg(x+w,y+h)),               # 1
      Part.LineSegment(self._from_svg(x+w,y+h), self._from_svg(x  ,y+h)),               # 2
      Part.LineSegment(self._from_svg(x  ,y+h), self._from_svg(x  ,y  )),               # 3
      # construction four corners
      Part.LineSegment(self._from_svg(x+rx,y   ), self._from_svg(x+rx,y+ry)),           # 4
      Part.LineSegment(self._from_svg(x+rx,y+ry), self._from_svg(x   ,y+ry)),           # 5
      Part.LineSegment(self._from_svg(x+w-rx,y   ), self._from_svg(x+w-rx,y+ry)),       # 6
      Part.LineSegment(self._from_svg(x+w-rx,y+ry), self._from_svg(x+w   ,y+ry)),       # 7
      Part.LineSegment(self._from_svg(x+w-rx,y+h   ), self._from_svg(x+w-rx,y+h-ry)),   # 8
      Part.LineSegment(self._from_svg(x+w-rx,y+h-ry), self._from_svg(x+w   ,y+h-ry)),   # 9
      Part.LineSegment(self._from_svg(x+rx,y+h   ), self._from_svg(x+rx,y+h-ry)),       # 10
      Part.LineSegment(self._from_svg(x+rx,y+h-ry), self._from_svg(x   ,y+h-ry))        # 11
    ], True)
    self.ske.addConstraint([
      ## outer construction corners
      Sketcher.Constraint('Coincident', i+0,2, i+1,1),
      Sketcher.Constraint('Coincident', i+1,2, i+2,1),
      Sketcher.Constraint('Coincident', i+2,2, i+3,1),
      Sketcher.Constraint('Coincident', i+3,2, i+0,1),
      ## inner construction corners
      Sketcher.Constraint('Coincident', i+4,2, i+5,1),
      Sketcher.Constraint('Coincident', i+6,2, i+7,1),
      Sketcher.Constraint('Coincident', i+8,2, i+9,1),
      Sketcher.Constraint('Coincident', i+10,2, i+11,1),
      ## inner construction equality
      Sketcher.Constraint('Equal', i+4, i+6),
      Sketcher.Constraint('Equal', i+4, i+8),
      Sketcher.Constraint('Equal', i+4, i+10),
      Sketcher.Constraint('Equal', i+5, i+7),
      Sketcher.Constraint('Equal', i+5, i+9),
      Sketcher.Constraint('Equal', i+5, i+11),
      ## corner cube outlines construction
      Sketcher.Constraint('PointOnObject',i+4,1, i+0),
      Sketcher.Constraint('PointOnObject',i+5,2, i+3),
      Sketcher.Constraint('PointOnObject',i+6,1, i+0),
      Sketcher.Constraint('PointOnObject',i+7,2, i+1),
      Sketcher.Constraint('PointOnObject',i+8,1, i+2),
      Sketcher.Constraint('PointOnObject',i+9,2, i+1),
      Sketcher.Constraint('PointOnObject',i+10,1, i+2),
      Sketcher.Constraint('PointOnObject',i+11,2, i+3),
      ## horizontal construction
      Sketcher.Constraint('Parallel', i, i+2),
      Sketcher.Constraint('Parallel', i, i+5),
      Sketcher.Constraint('Parallel', i, i+7),
      Sketcher.Constraint('Parallel', i, i+9),
      Sketcher.Constraint('Parallel', i, i+11),
      ## vertical construction
      Sketcher.Constraint('Parallel', i+1, i+3),
      Sketcher.Constraint('Parallel', i+1, i+4),
      Sketcher.Constraint('Parallel', i+1, i+6),
      Sketcher.Constraint('Parallel', i+1, i+8),
      Sketcher.Constraint('Parallel', i+1, i+10)
    ])
    ske.addGeometry([
      # sides of the rect
      Part.LineSegment(self._from_svg(x+rx  ,y     ), self._from_svg(x+w-rx,y     )),   # 12
      Part.LineSegment(self._from_svg(x+w   ,y+ry  ), self._from_svg(x+w   ,y+h-ry)),   # 13
      Part.LineSegment(self._from_svg(x+w-rx,y+h   ), self._from_svg(x+rx  ,y+h   )),   # 14
      Part.LineSegment(self._from_svg(x     ,y+h-ry), self._from_svg(x     ,y+ry  ))    # 15
    ])
    # arcs top left, top right, botton right, bottom left.
    # circles rotate counter clockwise. pi/2 is north, pi is west, 2*pi is east
    a_tl = self.objArc("", x+rx  , y+ry,   rx, ry, -2/2.*math.pi, -1/2.*math.pi, False, node, mat)
    a_tr = self.objArc("", x-rx+w, y+ry,   rx, ry, -1/2.*math.pi,  0/2.*math.pi, False, node, mat)
    a_br = self.objArc("", x-rx+w, y-ry+h, rx, ry,  0/2.*math.pi,  1/2.*math.pi, False, node, mat)
    a_bl = self.objArc("", x+rx  , y-ry+h, rx, ry,  1/2.*math.pi,  2/2.*math.pi, False, node, mat)
    if True:
      self.ske.addConstraint([
        # connect the corners to the edges. smooth
        Sketcher.Constraint('Tangent', a_tl[1][0],a_tl[1][1], i+12,1),
        Sketcher.Constraint('Tangent', i+12,2, a_tr[0][0],a_tr[0][1]),
        Sketcher.Constraint('Tangent', a_tr[1][0],a_tr[1][1], i+13,1),
        Sketcher.Constraint('Tangent', i+13,2, a_br[0][0],a_br[0][1]),
        Sketcher.Constraint('Tangent', a_br[1][0],a_br[1][1], i+14,1),
        Sketcher.Constraint('Tangent', i+14,2, a_bl[0][0],a_bl[0][1]),
        Sketcher.Constraint('Tangent', a_bl[1][0],a_bl[1][1], i+15,1),
        Sketcher.Constraint('Tangent', i+15,2, a_tl[0][0],a_bl[0][1]),
      ])
    if False:
      self.ske.addConstraint([
        # stitch the rounded rect to the construction grid
        Sketcher.Constraint('Coincident', i+12,1, i+4,1),
        Sketcher.Constraint('Coincident', i+12,2, i+6,1),
        Sketcher.Constraint('Coincident', i+13,1, i+7,2),
        Sketcher.Constraint('Coincident', i+13,2, i+9,2),
        Sketcher.Constraint('Coincident', i+14,1, i+8,1),
        Sketcher.Constraint('Coincident', i+14,2, i+10,1),
        Sketcher.Constraint('Coincident', i+15,1, i+11,2),
        Sketcher.Constraint('Coincident', i+15,2, i+5,2)
      ])
    if False and a_tr[3] is not None:         # ArcOfCirle has no majAxis
      self.ske.addConstraint([
        # make all major axis parallel, and same length
        Sketcher.Constraint('Parallel', a_tr[3], a_tl[3]),
        Sketcher.Constraint('Equal',    a_tr[3], a_tl[3]),
        Sketcher.Constraint('Parallel', a_tr[3], a_bl[3]),
        Sketcher.Constraint('Equal',    a_tr[3], a_bl[3]),
        Sketcher.Constraint('Parallel', a_tr[3], a_br[3]),
        # Sketcher.Constraint('Equal',    a_tr[3], a_br[3])     # makes everything immobole
      ])

    self.stats['objRect'] += 1


  def objRect(self, x, y, w, h, node, mat):
    self._matrix_from_svg(mat)
    i = int(self.ske.GeometryCount)
    self.ske.addGeometry([
      Part.LineSegment(self._from_svg(x  ,y  ), self._from_svg(x+w,y  )),
      Part.LineSegment(self._from_svg(x+w,y  ), self._from_svg(x+w,y+h)),
      Part.LineSegment(self._from_svg(x+w,y+h), self._from_svg(x  ,y+h)),
      Part.LineSegment(self._from_svg(x  ,y+h), self._from_svg(x  ,y  ))
    ], False)
    self.ske.addConstraint([
      Sketcher.Constraint('Coincident', i+0,2, i+1,1),
      Sketcher.Constraint('Coincident', i+1,2, i+2,1),
      Sketcher.Constraint('Coincident', i+2,2, i+3,1),
      Sketcher.Constraint('Coincident', i+3,2, i+0,1),
      Sketcher.Constraint('Parallel', i+2, i+0),
      Sketcher.Constraint('Parallel', i+3, i+1)
    ])
    self.stats['objRect'] += 1


  def objEllipse(self, cx, cy, rx, ry, node, mat):
    """
    We distinguish two cases. If it looks like a circle (after transformation),
    we produce a Circle, else we produce an Ellipse.
    The difference is clearly visible as we exposeInternalGeometry() of the Ellipse.
    """
    ### CAUTION: Keep in sync with objArc() below.
    self._matrix_from_svg(mat)
    c = self._from_svg(cx, cy, bbox=False)
    ori = self._from_svg(0, 0, bbox=False)
    vrx = self._from_svg(rx, 0, bbox=False) - ori
    vry = self._from_svg(0, ry, bbox=False) - ori
    i = int(self.ske.GeometryCount)

    if abs(vrx.Length - vry.Length) < epsilon:
      # it is a circle.
      self.bbox.add(c+vrx+vry)
      self.bbox.add(c-vrx-vry)
      self.ske.addGeometry([ Part.Circle(Center=c, Normal=Base.Vector(0,0,1), Radius=vrx.Length) ])
      self.stats['objEllipse'] += 1
    else:
      # major axis is defined by Center and S1,
      # major radius is the distance between Center and S1,
      # minor radius is the distance between S2 and the major axis.
      s1 = self._from_svg(cx+rx, cy, bbox=False)
      s2 = self._from_svg(cx, cy+ry, bbox=False)
      (V1,V2,V3,V4) = self._ellipse_vertices2d(c, s1, s2)
      self.bbox.add(V1)
      self.bbox.add(V2)
      self.bbox.add(V3)
      self.bbox.add(V4)
      self.ske.addGeometry([ Part.Ellipse(S1=V1, S2=V2, Center=c), ])
      if self.expose_int: self.ske.exposeInternalGeometry(i)
      self.stats['objEllipse'] += 1


  def objArc(self, d, cx, cy, rx, ry, st, en, closed, node, mat):
    """
    We ignore the path d, and produce a nice arc object.
    We distinguish two cases. If it looks like a circle, we produce ArcOfCircle,
    else we produce ArcOfEllipse.

    To find the arc end points, we use the value() property of the Circle and Ellipse
    objects. With Circle we have to take care of mirror and rotation ourselves.

    If closed, we connect the two ends to the center with lines. The connections
    are secured with constraints.
    Radii are currently not secured with constraints.
    """
    ### CAUTION: Keep in sync with objEllipse() above.
    # print("objArc: st,en,closed", st, en, closed)
    self._matrix_from_svg(mat)
    c = self._from_svg(cx, cy, bbox=False)
    ori = self._from_svg(0, 0, bbox=False)
    vrx = self._from_svg(rx, 0, bbox=False) - ori
    vry = self._from_svg(0, ry, bbox=False) - ori
    i = self.ske.GeometryCount
    (st_idx,en_idx) = (1,2)
    majAxisIdx = None

    if abs(vrx.Length - vry.Length) < epsilon:
      # it is a circle.
      self.bbox.add(c+vrx+vry)
      self.bbox.add(c-vrx-vry)
      ce = Part.Circle(Center=c, Normal=Base.Vector(0,0,1), Radius=vrx.Length)

      # Circles are immune to rotation and mirorring. Apply this manually.
      if self.yflip:
        (st,en) = (-en,-st)                     ## coord system is mirrored.
      else:
        (st,en) = (en,st)                       ## hmm.
        print("FIXME: ArcOfCircle() with yflip=False needs debugging.")
      r = Base.Matrix()
      r.rotateZ(st)
      pst = r.multiply(vrx)
      st = pst.getAngle(Base.Vector(1,0,0))     # ce.rotateZ() is a strange beast.
      pst = pst + c
      r = Base.Matrix()
      r.rotateZ(en)
      pen = r.multiply(vrx)
      en = pen.getAngle(Base.Vector(1,0,0))     # ce.rotateZ() is a strange beast.
      pen = pen + c

      self.ske.addGeometry([ Part.ArcOfCircle(ce, st, en) ])
      self.stats['objArc'] += 1

    else:
      # major axis is defined by Center and S1,
      # major radius is the distance between Center and S1,
      # minor radius is the distance between S2 and the major axis.
      s1 = self._from_svg(cx+rx, cy, bbox=False)
      s2 = self._from_svg(cx, cy+ry, bbox=False)
      (V1,V2,V3,V4) = self._ellipse_vertices2d(c, s1, s2)
      self.bbox.add(V1)
      self.bbox.add(V2)
      self.bbox.add(V3)
      self.bbox.add(V4)
      i = int(self.ske.GeometryCount)
      ce = Part.Ellipse(S1=V1, S2=V2, Center=c)
      self.ske.addGeometry([ Part.ArcOfEllipse(ce, st, en) ])
      if self.expose_int: self.ske.exposeInternalGeometry(i)
      majAxisIdx = i+1          # CAUTION: is that a safe assumption?
      self.stats['objArc'] += 1
      ## CAUTION: with yflip=True sketcher reverses the endpoints of
      ##          an ArcOfEllipse to: en=1, st=2
      ##          ArcOfCircle seems unaffected.
      if self.yflip: (st_idx,en_idx) = (2,1)
      r = Base.Matrix()
      r.rotateZ(st)
      pst = r.multiply(vrx) + c
      r = Base.Matrix()
      r.rotateZ(en)
      pen = r.multiply(vrx) + c

    j = self.ske.GeometryCount
    if closed:
      self.ske.addGeometry([
        Part.LineSegment(ce.value(en),c),
        Part.LineSegment(c,ce.value(st)) ])

      if True:          # when debugging deformations, switch off constriants first.
        self.ske.addConstraint([
          Sketcher.Constraint('Coincident', i+0,en_idx, j+0,1),   # arc with line
          Sketcher.Constraint('Coincident', j+1,2, i+0,st_idx),   # line with arc
          Sketcher.Constraint('Coincident', j+0,2, j+1,1),        # line with line
          Sketcher.Constraint('Coincident', j+0,2, i+0,3) ])      # line with center

    if False:    # some debugging circles.
      self.ske.addGeometry([
        # Part.Circle(Center=pst, Normal=Base.Vector(0,0,1), Radius=2),
        # Part.Circle(Center=pen, Normal=Base.Vector(0,0,1), Radius=3),
        # Part.Circle(Center=ce.value(st), Normal=Base.Vector(0,0,1), Radius=4),
        Part.Circle(Center=ce.value(en), Normal=Base.Vector(0,0,1), Radius=5)
        ], True)

    # we return the start, end and center points, as triple (sketcher_index, sketcher_index_point, Vector)
    return ((i+0, st_idx, ce.value(st)), (i+0, en_idx, ce.value(en)), (i+0, 3, c), majAxisIdx)


fcdoc = FreeCAD.newDocument(docname)
ske = fcdoc.addObject('Sketcher::SketchObject', 'Sketch_'+docname)

svg = InkSvg(pathgen=SketchPathGen(ske, yflip=True, expose_int=options.expose_internal_geometry))
svg.load(svgfile)       # respin of inkex.affect()
svg.traverse(options.ids)

if verbose > 0:
  print("svg2fcsketch %s, InkSvg %s" % (__version__, InkSvg.__version__))
if verbose >= 0:
  print(svg.stats)
# print(svg.docTransform, svg.docWidth, svg.docHeight, svg.dpi)

# print(ske)
#fcdoc.recompute()

if not options.outfile:
  import tempfile
  fcstdfile = tempfile.mktemp(prefix=docname, suffix='.fcstd')

fcdoc.saveAs(fcstdfile)
## Add GuiDocument.xml to the zip archive of fcstdfile
## to switch on default visibilitiy, and set a default camera.
camera_xml = ''
if True:        # switch off, if this causes errors. Nice to have.
  bb = svg.pathgen.bbox
  if verbose >= 0: print(bb)
  cx = bb.Center.x   # 35.246845 # bbox center
  cy = bb.Center.y   # 37.463238 # bbox center
  cz = bb.DiagonalLength * 0.5  # 51.437702 # focal distance: 1/2 of bbox diagonal
  zd = cz * 0.001             # 0.05 # +/- for far/near distance
  camera_xml = """<Camera settings="  OrthographicCamera { viewportMapping ADJUST_CAMERA position %f %f %f orientation 0 0 1  0 nearDistance %f farDistance %f aspectRatio 1 focalDistance %f height %f } "/>""" % (cx,cy,cz, cz-zd, cz+zd, cz, 2*cz)
guidoc_xml = """<?xml version='1.0' encoding='utf-8'?>
<Document SchemaVersion="1"><!-- as seen in FreeCAD 0.17 -->
    <ViewProviderData Count="1">
        <ViewProvider name="%s" expanded="0">
            <Properties Count="1">
                <Property name="Visibility" type="App::PropertyBool">
                    <Bool value="true"/>
                </Property>
            </Properties>
        </ViewProvider>
    </ViewProviderData>
    %s
</Document>
""" % ('Sketch_'+docname, camera_xml)

try:
  import zipfile
  z = zipfile.ZipFile(fcstdfile, 'a')
  z.writestr('GuiDocument.xml', guidoc_xml)
  z.close()
except:
  print(guidoc_xml)
  print("Warning: Failed to add GuiDocument.xml to %s -- camera and visibility are undefined." % fcstdfile)

if verbose > -1:
  print("%s written." % fcstdfile)

if not options.outfile:
  sys.stdout.write(open(fcstdfile).read())
  os.unlink(fcstdfile)
