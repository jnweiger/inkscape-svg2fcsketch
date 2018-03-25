#! /usr/bin/python
#
# (C) 2018 juergen@fabmail.org
# Distribute under GPL-2.0 or ask.
#
# References:
#  https://github.com/FreeCAD/FreeCAD/blob/master/src/Mod/Sketcher/TestSketcherApp.py
#  /usr/lib/freecad-daily/Mod/Sketcher/SketcherExample.py
#  https://en.wikipedia.org/wiki/Rytz%27s_construction#Computer_aided_solution
#
# v0.1 jw, initial draft refactoring inksvg to make it fit here.
# v0.2 jw, Introducing class SketchPathGen to seperate the sketch generator from the svg parser.
# v0.3 jw, correct _coord_from_svg() size and offset handling. Suppress
#          silly version printing, that would ruin an inkscape extension.
# V0.4 jw, Added GuiDocument.xml for visibility and camera defaults.
#          Using BoundBox() to compute camera placement.
# V0.5 jw, objEllipse() done correctly with _ellipse_vertices2d()
# V0.6 jw, objArc() done. ArcOfCircle() is a strange beast with rotation and mirroring.
#

from optparse import OptionParser
import os, sys, math, re

sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')

verbose=0       # 0=quiet, 1=normal
epsilon = 0.0000001

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

__version__ = '0.6'

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
parser.add_option('--smoothness', dest='smoothness',
                        type='float', default=float(0.2), action='store',
                        help='Curve smoothing (less for more [0.0001 .. 5]). Default: 0.2')
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


  def pathString(self, d, node, mat):
    """
    d is expected formatted as an svg path string here.
    """
    print("not impl. pathString: ", d, node, mat)


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
    print("not impl. objRoundedRect: ", x, y, w, h, rx, ry, node, mat)


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

