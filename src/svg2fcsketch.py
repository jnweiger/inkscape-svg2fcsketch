#! /usr/bin/python
#
# (C) 2018 juergen@fabmail.org
# Distribute under GPL-2.0 or ask.
#
# FROM:
#  https://github.com/FreeCAD/FreeCAD/blob/master/src/Mod/Sketcher/TestSketcherApp.py
#  /usr/lib/freecad-daily/Mod/Sketcher/SketcherExample.py
#
# v0.1 jw, initial draft refactoring inksvg to make it fit here.
# v0.2 jw, Introducing class SketchPathGen to seperate the sketch generator from the svg parser.

import os, sys, math, re
from optparse import OptionParser


sys.path.append('/usr/lib/freecad-daily/lib/')  # prefer daily over normal.
sys.path.append('/usr/lib/freecad/lib/')
from FreeCAD import Base
import Part, Sketcher
import ProfileLib.RegularPolygon as Poly

## INLINE_BLOCK_START
# for easier distribution, our Makefile can inline these imports
# sys.path.append('/home/src/github/jnweiger/inkscape-thunderlaser/src/')
sys.path.append(os.path.abspath(os.path.dirname(__file__))+'/../inksvg/src/')
from inksvg import InkSvg, PathGenerator
## INLINE_BLOCK_END

__version__ = '0.2'

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
  def __init__(self, ske, scale=1.0, yflip=True, yoff=0.0, stats=None):
    self.ske = ske
    self.m = None
    self.scale = scale
    self.yflip = yflip
    self.yoff = yoff

    # count how many objects of each type we generate.
    self.stats = {}
    if stats is not None: self.stats = stats
    for s in ('pathList', 'pathString', 'objRect', 'objRoundedRect',
              'objArc', 'objEllipse'):
      self.stats[s] = 0

  def _coord_from_svg(self, m=None):
    xs = self.scale
    ys = self.scale
    if self.yflip: ys = -self.scale
    if m is None: m = self.m
    m.move(0, self.yoff, 0)
    m.scale(xs, ys, 1)
    return m

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
    return m.multiply(Base.Vector(x, y, 0))


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
    i = int(ske.GeometryCount)      # 0
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
    ske.addGeometry([
      Part.LineSegment(self._from_svg(x  ,y  ), self._from_svg(x+w,y  )),
      Part.LineSegment(self._from_svg(x+w,y  ), self._from_svg(x+w,y+h)),
      Part.LineSegment(self._from_svg(x+w,y+h), self._from_svg(x  ,y+h)),
      Part.LineSegment(self._from_svg(x  ,y+h), self._from_svg(x  ,y  ))
    ], False)
    ske.addConstraint([
      Sketcher.Constraint('Coincident',i+0,2, i+1,1),
      Sketcher.Constraint('Coincident',i+1,2, i+2,1),
      Sketcher.Constraint('Coincident',i+2,2, i+3,1),
      Sketcher.Constraint('Coincident',i+3,2, i+0,1),
    ])
    self.stats['objRect'] += 1

  def objEllipse(self, cx, xy, rx, ry, node, mat):
    print("not impl. objEllipse: ", cx, xy, rx, ry, node, mat)

  def objArc(self, d, cx, cy, rx, ry, st, en, cl, node, mat):
    """
    We ignore the path d, and produce a nice arc object.
    """


fcdoc = FreeCAD.newDocument(docname)
ske = fcdoc.addObject('Sketcher::SketchObject', 'Sketch_'+docname)

gen_stats = {}
svg = InkSvg(pathgen=SketchPathGen(ske, stats=gen_stats, scale=1.0, yflip=True, yoff=-100.0))
svg.load(svgfile)       # respin of inkex.affect()
selected = svg.getElementsByIds(options.ids)

## FIXME: hide this in svg.traverse([ids...])
if len(selected):
  # Traverse the selected objects
  for node in selected:
    transform = svg.recursivelyGetEnclosingTransform(node)
    svg.recursivelyTraverseSvg([node], transform)
else:
  # Traverse the entire document building new, transformed paths
  svg.recursivelyTraverseSvg(svg.document.getroot(), svg.docTransform)

print(ske)
print(gen_stats)

#fcdoc.recompute()

fcdoc.saveAs(fcstdfile)
print("%s written." % fcstdfile)

