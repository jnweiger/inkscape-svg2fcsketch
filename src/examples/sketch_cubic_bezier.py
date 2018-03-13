#!/usr/bin/python
# References: 
# - https://en.wikipedia.org/wiki/Composite_B%C3%A9zier_curve
# - https://en.wikipedia.org/wiki/B-spline#Relationship_to_piecewise/composite_B%C3%A9zier
#
# A series of SVG curveto (C) path commands specifies a composite bezier spline of cubic or 3rd order bezier curves.
#       The endpoints of each curveto can have C0, C1, or C2 continuity, represented in inkscape by handles of
#       arbitrary length and angle (C0), handles with a 180° restriction (C1), or handles with 180° and same length restriction (C2).
#       An addintional path attribute sodipodi:nodetypes="ccszzcsc" is used to specify C0=c, C1=s or C2=z restrictions. 
#       There is one letter per node. Nodes of polygon segments (L, H, V) have letter c.
#
# Freecad sketches support higher order B-Splines. We can construct an exact representation of SVG curveto paths with C0, C1, or C2 
#       properties by composing cubic bezier splines. The handles as known from inkscape are represented by lines in construction mode
#       connecting curve points (two end points per cubic bezier spline) with their control points (two inner points per cubic bezier spline).
#       
# An inkscape elliptic arc path command A has 4 groups of x,y 0 0 1 x,y coordinates, a total of 28. 
#       TODO: understand why 4 groups, conversion to one of the freecad arc types. 
#             Probably easier done by evaluating the sodipodi:type,cx,xy,rx,ry,start,end,open
#             In any case, transform attributes must be calculated, as freecad sketches have no transforms for their elements.
#             Can inksvg.py resolve all that into proper splines?
 
import os, sys, math
sys.path.append('/usr/share/inkscape/extensions')       # Linux
import inkex

# Inkscape uses a well defined subset of the svg specification.
# It is probably easier done as an inkscape exporter than as a freecad importer.

# fcstd files for FreeCAD are zip archives containing at least a Document.xml file. 
# An additional GuiDocument.xml file makes the sketch visible when loading. 
# - Otherwise the 'center all' tool (V,F) can be used to make things visible.

