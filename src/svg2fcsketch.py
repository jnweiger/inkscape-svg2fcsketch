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

verbose=0       # 0=quiet, 1=normal
epsilon = 0.00001

if verbose == 0:
  os.dup2(1,99)         # hack to avoid silly version string printing.
  f = open("/dev/null", "w")
  os.dup2(f.fileno(), 1)

# The for version printing code has
# src/App/Application.cpp: if (!(mConfig["Verbose"] == "Strict"))
# but we cannot call SetConfig('Verbose', 'Strict') early enough.
from FreeCAD import Base, BoundBox
sys.stdout.flush()      # push silly version string into /dev/null

if verbose == 0:
  f.close()
  os.dup2(99,1)         # back in cansas.

import Part, Sketcher   # causes SEGV if Base is not yet imported from FreeCAD
import ProfileLib.RegularPolygon as Poly

## INLINE_BLOCK_START
# for easier distribution, our Makefile can inline these imports
# sys.path.append('/home/src/github/jnweiger/inkscape-thunderlaser/src/')
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/../inksvg/src/')
from inksvg import InkSvg, PathGenerator
## INLINE_BLOCK_END

__version__ = '0.7'

parser = OptionParser(usage="\n    %prog [options] SVGFILE [OUTFILE]\n\nTry --help for details.")
parser.add_option("-o", "--outfile", dest="outfile",
                         help="write fcstd to OUTPUT. Default: same as input file, but with .fcstd suffix instead.", metavar="OUTPUT")
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


class SketchPathGen(PathGenerator):
  """
  Generate XML code for a FreeCAD sketch.
  """
  def __init__(self, ske, yflip=True):
    self.ske = ske
    self.m = None
    self.yflip = yflip
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

    """
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
    for sp in path:
      # These are the off by one's: four points -> three lines -> two constraints.
      prev_idx = None
      j = 0
      while j < len(sp)-1:
        # lists of three, http://wiki.inkscape.org/wiki/index.php/Python_modules_for_extensions#cubicsuperpath.py
        (h0,p1,h1) = sp[j]      
        j = j+1
        while j < len(sp):
          (h2,p2,h3) = sp[j]
          if not self._same_point(p1, p2):
            break               # no null-segments, please!
          j += 1
        if j >= len(sp):
          break                 # nothing left.
        
        if self._same_point(h1, p1) and self._same_point(h2, p2):
          # it is a straigth line
          i = int(self.ske.GeometryCount)
          self.ske.addGeometry([Part.LineSegment(self._from_svg(p1[0], p1[1]), self._from_svg(p2[0], p2[1]))])
        else:
          # it is a spline
          print("spline half. control points not visualized, not constaint", p1, h1, h2, p2)
          i = int(self.ske.GeometryCount)
          self.ske.addGeometry([Part.BSplineCurve([
            self._from_svg(p1[0], p1[1]),
            self._from_svg(h1[0], h1[1]),
            self._from_svg(h2[0], h2[1]),
            self._from_svg(p2[0], p2[1])
            ], None,None,False,3,None,False),False])
          self.ske.exposeInternalGeometry(i)
        # more internal geometry ...
        if prev_idx is not None:
          self.ske.addConstraint([Sketcher.Constraint('Coincident', prev_idx,2, i,1)])
        prev_idx = i
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

    print("objRoundedRect: ", x, y, w, h, rx, ry, node.get('id'), mat)

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
    if False:
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
      self.ske.exposeInternalGeometry(i)
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
      self.ske.exposeInternalGeometry(i)
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

svg = InkSvg(pathgen=SketchPathGen(ske, yflip=True))
svg.load(svgfile)       # respin of inkex.affect()
svg.traverse(options.ids)

if verbose > 0:
  print("svg2fcsketch %s, InkSvg %s" % (__version__, InkSvg.__version__))
print(svg.stats)
# print(svg.docTransform, svg.docWidth, svg.docHeight, svg.dpi)

# print(ske)
#fcdoc.recompute()

fcdoc.saveAs(fcstdfile)
## Add GuiDocument.xml to the zip archive of fcstdfile
## to switch on default visibilitiy, and set a default camera.
camera_xml = ''
if True:        # switch off, if this causes errors. Nice to have.
  bb = svg.pathgen.bbox
  print(bb)
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

