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
  def __init__(self, ske):
    self.ske = ske

  def pathVertices(self, d, node, mat):
    """
    d is expected formatted as an svg path string here.
    """
    print("not impl. pathVertices: ", d, node, mat)

  def simplePath(self, d, node, mat):
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

  def roundedRect(self, x, y, w, h, rx, ry, node, mat):
    print("not impl. roundedRect: ", x, y, w, h, rx, ry, node, mat)



fcdoc = FreeCAD.newDocument(docname)
ske = fcdoc.addObject('Sketcher::SketchObject', 'Sketch_'+docname)

svg = InkSvg(pathgen=SketchPathGen(ske))
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

print(svg.paths)        # do not use svg.paths here. Use the ske object!

##ske.addConstraint(con)
#doc.recompute()

#doc.saveAs(fcstdfile)
#print("%s written." % fcstdfile)

