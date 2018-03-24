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

__version__ = '0.4'

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
    seen in https://en.wikipedia.org/wiki/Rytz%27s_construction
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


  def _from_svg(self, x, y, m=None):
    """
    Converts a 2D Vector from SVG into a FreeCAD 3D vector applying Base.Matrix m.
    """
    if m is None: m = self.m
    v = m.multiply(Base.Vector(x, y, 0))
    self.bbox.add(v)
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
    self._matrix_from_svg(mat)
    print("svg mat", mat)
    print('elli2 transform="matrix(1,0,-0.26794919,1,79.58091,0)">')
    print("fc mat", self.m)
    c = self._from_svg(cx, cy)
    ori = self._from_svg(0, 0)
    vrx = self._from_svg(rx, 0) - ori
    vry = self._from_svg(0, ry) - ori
    if abs(vrx.Length - vry.Length) < epsilon:
      # it is a circle.
      self.ske.addGeometry([ Part.Circle(Center=c, Normal=Base.Vector(0,0,1), Radius=vrx.Length) ])
      self.stats['objEllipse'] += 1
    else:
      if True:
        if rx > ry:
          e = Part.Ellipse(c, rx, ry)
        else:
          e = Part.Ellipse(c, ry, rx)
        (tr,sc,ro) = self._decompose_matrix2d()
        m = Base.Matrix()
        m.rotateZ(ro[0])
        e.rotate(Base.Placement(m))
        self.ske.addGeometry([ e ])
        self.ske.exposeInternalGeometry(self.ske.GeometryCount-1)
        s1 = self._from_svg(cx+rx, cy)
        s2 = self._from_svg(cx, cy+ry)
        self.ske.addGeometry([ Part.Circle(Center=s1, Normal=Base.Vector(0,0,1), Radius=4),
                               Part.Circle(Center=s2, Normal=Base.Vector(0,0,1), Radius=6)
           ], True)
        self.stats['objEllipse'] += 1
      else:
        # major axis is defined by Center and S1,
        # major radius is the distance between Center and S1,
        # minor radius is the distance between S2 and the major axis.
        s1 = self._from_svg(cx+rx, cy)
        s2 = self._from_svg(cx, cy+ry)
        print("s1 ", s1)
        print("s2 ", s2)
        print("c  ", c)
        (tr,sc,ro) = self._decompose_matrix2d()
        print("tr ", tr)
        print("sc ", sc)
        print("ro ", ro)
        print("FIXME: Wrong computation, the ellipsis is not long enough, and it is too fat.")
        i = int(self.ske.GeometryCount)
        self.ske.addGeometry([ Part.Ellipse(S1=s1, S2=s2, Center=c), ])
        self.ske.exposeInternalGeometry(self.ske.GeometryCount-1)
        self.ske.addGeometry([ Part.Circle(Center=s1, Normal=Base.Vector(0,0,1), Radius=4),
                               Part.Circle(Center=s2, Normal=Base.Vector(0,0,1), Radius=6)
           ], True)
        self.stats['objEllipse'] += 1


  def objArc(self, d, cx, cy, rx, ry, st, en, cl, node, mat):
    """
    We ignore the path d, and produce a nice arc object.
    """


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

