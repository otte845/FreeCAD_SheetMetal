########################################################################
#
#  SheetMetalLoft.py
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2 of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
########################################################################
"""Sheet metal loft between two wires (open or closed).

The two profiles are converted to polylines and the space between them is
filled with planar triangles. Closed profiles are cut open at a start
vertex and then treated like open ones.
"""

import math
import os

import FreeCAD
import Part

import SheetMetalTools

translate = FreeCAD.Qt.translate
icons_path = SheetMetalTools.icons_path

# Geometric tolerance used to detect coincident points / degenerate triangles.
LOFT_TOLERANCE = 1e-6
# A candidate triangle whose normal differs from the previous triangle's
# normal by more than acos(FOLD_BACK_COS) is considered folded back.
FOLD_BACK_COS = -0.5


###################################################################################################
# Profile preparation
###################################################################################################

def getProfileWire(obj):
    """Return the single (open or closed) wire contained in the shape of `obj`."""
    shape = getattr(obj, "Shape", None)
    if shape is None or shape.isNull() or not shape.Edges:
        raise ValueError(f"'{obj.Label}' has no edges")
    groups = Part.sortEdges(shape.Edges)
    if len(groups) != 1:
        raise ValueError(f"'{obj.Label}' must contain exactly one connected wire "
                         f"(found {len(groups)})")
    return Part.Wire(groups[0])


def _isLineEdge(edge):
    return isinstance(edge.Curve, (Part.Line, Part.LineSegment))


def _edgeParameters(edge, start, numSegments):
    """Equal-length parameters along `edge`, ordered from vertex `start`."""
    forward = ((edge.valueAt(edge.FirstParameter) - start).Length
               <= (edge.valueAt(edge.LastParameter) - start).Length)
    if numSegments == 1:
        params = [edge.FirstParameter, edge.LastParameter]
    else:
        length = edge.Length
        params = [edge.FirstParameter]
        params += [edge.getParameterByLength(length * k / numSegments)
                   for k in range(1, numSegments)]
        params.append(edge.LastParameter)
    return params if forward else params[::-1]


def discretizeWire(wire, minLineLength):
    """Convert a wire to an ordered list of polyline points.

    Straight edges are kept as a single segment. Any other curve (arcs,
    splines, ...) is split into equal-length chords not shorter than
    `minLineLength` (at least one chord per edge, three for a closed edge
    such as a full circle). For a closed wire the last point repeats the
    first one.
    """
    edges = wire.OrderedEdges
    verts = [v.Point for v in wire.OrderedVertexes]
    if len(verts) == len(edges):
        verts.append(verts[0])
    points = [verts[0]]
    for k, edge in enumerate(edges):
        start, end = verts[k], verts[k + 1]
        if _isLineEdge(edge) or minLineLength <= 0:
            numSegments = 1
        else:
            numSegments = max(1, int(math.floor(edge.Length / minLineLength)))
        if edge.isClosed():
            numSegments = max(3, numSegments)
        params = _edgeParameters(edge, start, numSegments)
        for p in [edge.valueAt(q) for q in params[1:-1]] + [end]:
            if (p - points[-1]).Length > LOFT_TOLERANCE:
                points.append(p)
    if wire.isClosed():
        points[-1] = FreeCAD.Vector(points[0])
    return points


def isClosedPolyline(points):
    return len(points) > 3 and (points[-1] - points[0]).Length < LOFT_TOLERANCE


def polylineCenter(points):
    """Length weighted geometric center of a polyline."""
    total = 0.0
    center = FreeCAD.Vector()
    for a, b in zip(points[:-1], points[1:]):
        length = (b - a).Length
        center += (a + b) * (0.5 * length)
        total += length
    if total < LOFT_TOLERANCE:
        return FreeCAD.Vector(points[0])
    return center * (1.0 / total)


def windingVector(points, center):
    """Sum of (start - center) x segment over all polyline segments."""
    winding = FreeCAD.Vector()
    for a, b in zip(points[:-1], points[1:]):
        winding += (a - center).cross(b - a)
    return winding


def _windingSign(points, center, direction):
    """+1 / -1 if the winding agrees / disagrees with `direction`,
    0 if it can't be decided (straight wire, or winding normal to direction).
    """
    winding = windingVector(points, center)
    if winding.Length < LOFT_TOLERANCE or direction.Length < LOFT_TOLERANCE:
        return 0
    cosine = winding.dot(direction) / (winding.Length * direction.Length)
    if abs(cosine) < 1e-3:
        return 0
    return 1 if cosine > 0 else -1


def _endpointsMismatch(pts1, pts2):
    return (pts1[0] - pts2[0]).Length + (pts1[-1] - pts2[-1]).Length


def _orientationFlags(pts1, pts2):
    """Return (reverse1, reverse2) so both polylines run 'clockwise'."""
    c1, c2 = polylineCenter(pts1), polylineCenter(pts2)
    direction = c2 - c1
    sign1 = _windingSign(pts1, c1, direction)
    sign2 = _windingSign(pts2, c2, direction)
    rev1, rev2 = sign1 < 0, sign2 < 0
    o1 = pts1[::-1] if rev1 else pts1
    o2 = pts2[::-1] if rev2 else pts2
    if sign2 == 0:
        rev2 = _endpointsMismatch(o1, o2[::-1]) < _endpointsMismatch(o1, o2)
    elif sign1 == 0:
        rev1 = _endpointsMismatch(o1[::-1], o2) < _endpointsMismatch(o1, o2)
    return rev1, rev2


def orientProfiles(pts1, pts2):
    """Orient both polylines so they run 'clockwise' around the loft axis.

    The winding vector of each wire must point in the same direction as
    the vector from the center of wire 1 to the center of wire 2. When the
    winding can't decide (straight wires, coplanar profiles...) the wire
    is aligned to the other one by matching its end points.
    Returns new lists; the inputs are not modified.
    """
    rev1, rev2 = _orientationFlags(pts1, pts2)
    return (pts1[::-1] if rev1 else list(pts1)), (pts2[::-1] if rev2 else list(pts2))


DEFAULT_START_DIRECTION = FreeCAD.Vector(1, 0, 0)


def _profileDirection(direction, axis):
    """`direction` projected on the plane normal to the loft `axis`. Falls
    back to the global axes when `direction` is (nearly) along the axis."""
    if axis.Length > LOFT_TOLERANCE:
        axis = FreeCAD.Vector(axis).normalize()
    else:
        axis = None
    candidates = [direction, FreeCAD.Vector(1, 0, 0), FreeCAD.Vector(0, 1, 0),
                  FreeCAD.Vector(0, 0, 1)]
    for d in candidates:
        if d is None or d.Length < LOFT_TOLERANCE:
            continue
        d = FreeCAD.Vector(d).normalize()
        if axis is not None:
            d = d - axis * d.dot(axis)
        if d.Length > 1e-3:
            return d.normalize()
    return FreeCAD.Vector(1, 0, 0)


def startVertexByDirection(points, direction):
    """Index of the vertex of a closed polyline closest to the ray leaving
    the polyline center along `direction` (smallest angle to it)."""
    center = polylineCenter(points)
    best, bestCos = 0, None
    for k, p in enumerate(points[:-1]):
        v = p - center
        if v.Length < LOFT_TOLERANCE:
            continue
        cos = v.dot(direction) / v.Length
        if bestCos is None or cos > bestCos + 1e-9:
            best, bestCos = k, cos
    return best


def nearestVertex(points, target):
    """Index of the vertex of a closed polyline nearest to `target`."""
    loop = points[:-1]
    return min(range(len(loop)), key=lambda k: (loop[k] - target).Length)


def rotateClosedPolyline(points, start):
    """Closed polyline starting (and ending) at vertex `start`."""
    loop = points[:-1]
    loop = loop[start:] + loop[:start]
    return loop + [FreeCAD.Vector(loop[0])]


def alignStarts(pts1, pts2, direction=None):
    """Cut closed polylines open at matching start vertices.

    Both closed: wire 1 starts at its vertex closest to `direction` (seen
    from its center, projected normal to the loft axis), wire 2 at its
    vertex nearest to that one. One closed: it starts at its vertex nearest
    to the start of the open one. The polylines must already be oriented.
    """
    closed1, closed2 = isClosedPolyline(pts1), isClosedPolyline(pts2)
    if closed1:
        if closed2:
            axis = polylineCenter(pts2) - polylineCenter(pts1)
            if direction is None:
                direction = DEFAULT_START_DIRECTION
            start = startVertexByDirection(pts1, _profileDirection(direction, axis))
        else:
            start = nearestVertex(pts1, pts2[0])
        pts1 = rotateClosedPolyline(pts1, start)
    if closed2:
        pts2 = rotateClosedPolyline(pts2, nearestVertex(pts2, pts1[0]))
    return pts1, pts2


###################################################################################################
# Triangulation
###################################################################################################

def _cumulativeLengths(points):
    lengths = [0.0]
    for a, b in zip(points[:-1], points[1:]):
        lengths.append(lengths[-1] + (b - a).Length)
    total = lengths[-1]
    if total < LOFT_TOLERANCE:
        return [0.0] * len(points)
    return [s / total for s in lengths]


def triangleNormal(tri):
    a, b, c = tri
    return (b - a).cross(c - a)


def isDegenerateTriangle(tri, tol=LOFT_TOLERANCE):
    a, b, c = tri
    longest = max((b - a).Length, (c - b).Length, (a - c).Length)
    if longest < tol:
        return True
    # |normal| / longest edge == height over the longest edge.
    return triangleNormal(tri).Length / longest < tol


def isValidTriangle(tri, prevNormal, tol=LOFT_TOLERANCE):
    """A triangle is valid if it is not degenerate and does not fold back
    over the previously created triangle."""
    if isDegenerateTriangle(tri, tol):
        return False
    if prevNormal is not None:
        n = triangleNormal(tri)
        n.normalize()
        if n.dot(prevNormal) < FOLD_BACK_COS:
            return False
    return True


class LoftFront:
    """State of the advancing front while walking both polylines."""

    def __init__(self, pts1, pts2):
        self.pts1 = pts1
        self.pts2 = pts2
        self.s1 = _cumulativeLengths(pts1)
        self.s2 = _cumulativeLengths(pts2)
        self.i = 0
        self.j = 0
        self.prevNormal = None

    def candidate1(self, i=None, j=None):
        """Triangle advancing on wire 1: (P1[i], P1[i+1], P2[j])."""
        i = self.i if i is None else i
        j = self.j if j is None else j
        if i >= len(self.pts1) - 1:
            return None
        return (self.pts1[i], self.pts1[i + 1], self.pts2[j])

    def candidate2(self, i=None, j=None):
        """Triangle advancing on wire 2: (P1[i], P2[j+1], P2[j])."""
        i = self.i if i is None else i
        j = self.j if j is None else j
        if j >= len(self.pts2) - 1:
            return None
        return (self.pts1[i], self.pts2[j + 1], self.pts2[j])

    def finished(self):
        return self.i >= len(self.pts1) - 1 and self.j >= len(self.pts2) - 1


# Decision functions: called when both candidate triangles are equally valid.
# Signature: f(front, tri1, tri2) -> 1 to advance on wire 1, 2 for wire 2.

def chooseParametric(front, tri1, tri2):
    """Advance on the wire whose next point is less far along its own
    normalized length, so both wires progress at the same rate."""
    return 1 if front.s1[front.i + 1] <= front.s2[front.j + 1] else 2


def chooseShortestDiagonal(front, tri1, tri2):
    """Classic heuristic: keep the shortest new rung between the wires."""
    rung1 = (front.pts1[front.i + 1] - front.pts2[front.j]).Length
    rung2 = (front.pts2[front.j + 1] - front.pts1[front.i]).Length
    return 1 if rung1 <= rung2 else 2


def _unitNormal(tri):
    n = triangleNormal(tri)
    n.normalize()
    return n


def _bestFirstStep(front, tri1, tri2):
    """Look two steps ahead from the start: evaluate the four paths
    (first triangle x next triangle) and return the first step of the
    path with the smallest bend between its two triangles."""
    best, bestDot = None, None
    for pick, first in ((1, tri1), (2, tri2)):
        if isDegenerateTriangle(first):
            continue
        n1 = _unitNormal(first)
        i = front.i + (1 if pick == 1 else 0)
        j = front.j + (1 if pick == 2 else 0)
        for second in (front.candidate1(i, j), front.candidate2(i, j)):
            if second is None or not isValidTriangle(second, n1):
                continue
            dot = n1.dot(_unitNormal(second))
            if bestDot is None or dot > bestDot:
                best, bestDot = pick, dot
    return best


def chooseMinDihedral(front, tri1, tri2):
    """Choose the triangle that bends least relative to the previous one
    (smoothest / least concave surface)."""
    if front.prevNormal is None:
        pick = _bestFirstStep(front, tri1, tri2)
        return pick if pick is not None else chooseShortestDiagonal(front, tri1, tri2)
    n1 = _unitNormal(tri1)
    n2 = _unitNormal(tri2)
    return 1 if n1.dot(front.prevNormal) >= n2.dot(front.prevNormal) else 2


def chooseMaxMinAngle(front, tri1, tri2):
    """Choose the better shaped triangle (largest smallest inner angle)."""
    return 1 if _minAngle(tri1) >= _minAngle(tri2) else 2


def _minAngle(tri):
    angles = []
    for k in range(3):
        a, b, c = tri[k], tri[(k + 1) % 3], tri[(k + 2) % 3]
        u, v = b - a, c - a
        if u.Length < LOFT_TOLERANCE or v.Length < LOFT_TOLERANCE:
            return 0.0
        angles.append(u.getAngle(v))
    return min(angles)


TRIANGULATION_STRATEGIES = {
    "MinDihedral": chooseMinDihedral,
    "ShortestDiagonal": chooseShortestDiagonal,
    "Parametric": chooseParametric,
    "MaxMinAngle": chooseMaxMinAngle,
}
DEFAULT_STRATEGY = "MinDihedral"
# Closed profiles: largest allowed lead (fraction of the normalized length)
# of one wire over the other before the chooser is overruled. Prevents
# twisted lofts where one wire is walked around while the other waits.
CLOSED_MAX_LAG = 0.25


def _stepLag(front, pick):
    if pick == 1:
        return abs(front.s1[front.i + 1] - front.s2[front.j])
    return abs(front.s2[front.j + 1] - front.s1[front.i])


def triangulateSteps(pts1, pts2, chooser=chooseMinDihedral, maxLag=None):
    """Walk both polylines and return the non-degenerate triangles as
    (pick, i, j, triangle) tuples, where (i, j) is the front position
    before the step and pick tells which wire advanced (1 or 2)."""
    front = LoftFront(pts1, pts2)
    steps = []
    while not front.finished():
        tri1 = front.candidate1()
        tri2 = front.candidate2()
        if tri1 is not None and tri2 is not None:
            valid1 = isValidTriangle(tri1, front.prevNormal)
            valid2 = isValidTriangle(tri2, front.prevNormal)
            if valid1 != valid2:
                pick = 1 if valid1 else 2
            else:
                pick = chooser(front, tri1, tri2)
                other = 3 - pick
                if (maxLag is not None and _stepLag(front, pick) > maxLag
                        and _stepLag(front, other) < _stepLag(front, pick)):
                    pick = other
        else:
            pick = 1 if tri1 is not None else 2
        tri = tri1 if pick == 1 else tri2
        if not isDegenerateTriangle(tri):
            steps.append((pick, front.i, front.j, tri))
            front.prevNormal = _unitNormal(tri)
        if pick == 1:
            front.i += 1
        else:
            front.j += 1
    return steps


def triangulateProfiles(pts1, pts2, chooser=chooseMinDihedral):
    """Walk both polylines and return the list of triangles joining them."""
    bothClosed = isClosedPolyline(pts1) and isClosedPolyline(pts2)
    maxLag = CLOSED_MAX_LAG if bothClosed else None
    return [step[3] for step in triangulateSteps(pts1, pts2, chooser, maxLag)]


def makeTriangleShell(triangles):
    if not triangles:
        raise ValueError("Loft produced no triangles")
    faces = [Part.Face(Part.makePolygon([a, b, c, a])) for a, b, c in triangles]
    shell = Part.Shell(faces)
    shell.sewShape()
    if len(shell.Shells) == 1:
        return shell.Shells[0]
    return shell


def prepareProfiles(wire1, wire2, minLineLength, startDirection=None):
    """Discretize, orient and (for closed wires) align the start vertices."""
    pts1 = discretizeWire(wire1, minLineLength)
    pts2 = discretizeWire(wire2, minLineLength)
    pts1, pts2 = orientProfiles(pts1, pts2)
    return alignStarts(pts1, pts2, startDirection)


def smLoftShell(wire1, wire2, minLineLength=25.0, strategy=DEFAULT_STRATEGY,
                startDirection=None):
    """Create a triangulated shell between two wires."""
    pts1, pts2 = prepareProfiles(wire1, wire2, minLineLength, startDirection)
    chooser = TRIANGULATION_STRATEGIES.get(strategy, TRIANGULATION_STRATEGIES[DEFAULT_STRATEGY])
    return makeTriangleShell(triangulateProfiles(pts1, pts2, chooser))


###################################################################################################
# Feature object
###################################################################################################

class SMLoft:
    def __init__(self, obj, profiles):
        """Sheet metal loft between two wires (open or closed)."""
        self.addVerifyProperties(obj)
        obj.Profiles = profiles
        obj.Proxy = self

    def addVerifyProperties(self, obj):
        SheetMetalTools.smAddProperty(obj,
            "App::PropertyLinkList",
            "Profiles",
            translate("App::Property", "The two wire profiles to loft between"))
        SheetMetalTools.smAddLengthProperty(obj,
            "MinLineLength",
            translate("App::Property",
                      "Minimum length of the lines used to approximate curved segments"),
            25.0)
        SheetMetalTools.smAddEnumProperty(obj,
            "TriangulationMode",
            "Rule used to pick the next triangle when both candidates are valid",
            list(TRIANGULATION_STRATEGIES.keys()),
            DEFAULT_STRATEGY)
        SheetMetalTools.smAddProperty(obj,
            "App::PropertyVector",
            "StartDirection",
            translate("App::Property",
                      "Closed profiles: the first profile starts at its vertex closest to this "
                      "direction (seen from its center), the second at its vertex nearest "
                      "to that one"),
            FreeCAD.Vector(DEFAULT_START_DIRECTION))

    def execute(self, fp):
        self.addVerifyProperties(fp)
        if len(fp.Profiles) != 2:
            raise ValueError("SheetMetal Loft needs exactly two profiles")
        wire1 = getProfileWire(fp.Profiles[0])
        wire2 = getProfileWire(fp.Profiles[1])
        fp.Shape = smLoftShell(wire1, wire2,
                               minLineLength=fp.MinLineLength.Value,
                               strategy=fp.TriangulationMode,
                               startDirection=fp.StartDirection)


###################################################################################################
# Gui code
###################################################################################################

if SheetMetalTools.isGuiLoaded():
    Gui = FreeCAD.Gui


    class SMLoftViewProvider(SheetMetalTools.SMViewProvider):
        """Part / Part WB style ViewProvider."""

        def getIcon(self):
            return os.path.join(icons_path, "SheetMetal_Loft.svg")

        def claimChildren(self):
            if hasattr(self, "Object") and hasattr(self.Object, "Profiles"):
                return list(self.Object.Profiles)
            return []


    class AddLoftCommandClass:
        """Sheet metal loft command."""

        def GetResources(self):
            return {
                    "Pixmap": os.path.join(icons_path, "SheetMetal_Loft.svg"),
                    "MenuText": translate("SheetMetal", "Make Loft"),
                    "ToolTip": translate(
                        "SheetMetal",
                        "Create a sheet metal loft between two wires\n"
                        "1. Select two objects, each containing one open or closed wire.\n"
                        "2. Use Property editor to modify other parameters",
                        ),
                    }

        def Activated(self):
            profiles = Gui.Selection.getSelection()
            for profile in profiles:
                try:
                    getProfileWire(profile)
                except ValueError as e:
                    SheetMetalTools.SMLogger.error(str(e))
                    return
            doc = FreeCAD.ActiveDocument
            newObj, _ = SheetMetalTools.smCreateNewObject(profiles[0], "Loft",
                                                          allowPartDesign=False)
            if newObj is None:
                return
            SMLoft(newObj, profiles)
            SMLoftViewProvider(newObj.ViewObject)
            SheetMetalTools.smAddNewObject(profiles[0], newObj, None)
            profiles[1].ViewObject.Visibility = False
            doc.commitTransaction()

        def IsActive(self):
            sel = Gui.Selection.getSelection()
            return len(sel) == 2 and all(hasattr(o, "Shape") for o in sel)


    Gui.addCommand("SheetMetal_Loft", AddLoftCommandClass())
