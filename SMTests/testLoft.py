import math
import unittest

import FreeCAD
import Part

import SheetMetalLoft as Loft

V = FreeCAD.Vector


def _arcWire(radius, z, reverse=False):
    pts = [V(radius, 0, z), V(0, radius, z), V(-radius, 0, z)]
    if reverse:
        pts.reverse()
    return Part.Wire([Part.Arc(*pts).toShape()])


def _rectWire(dx, dy, z):
    return Part.makePolygon([V(dx, dy, z), V(-dx, dy, z), V(-dx, -dy, z), V(dx, -dy, z),
                             V(dx, dy, z)])


class TestSheetMetalLoft(unittest.TestCase):
    def test_lines_are_not_subdivided(self):
        wire = Part.makePolygon([V(0, 0, 0), V(100, 0, 0), V(100, 100, 0)])
        pts = Loft.discretizeWire(wire, 25.0)
        self.assertEqual(len(pts), 3)

    def test_arc_is_split_by_min_line_length(self):
        wire = _arcWire(50, 0)  # length = 50 * pi ~ 157
        pts = Loft.discretizeWire(wire, 25.0)
        self.assertEqual(len(pts), int(math.floor(50 * math.pi / 25.0)) + 1)
        for a, b in zip(pts[:-1], pts[1:]):
            self.assertGreaterEqual((b - a).Length, 25.0 * 0.95)

    def test_orientation_matches_reversed_profile(self):
        p1 = Loft.discretizeWire(_arcWire(50, 0), 25.0)
        p2 = Loft.discretizeWire(_arcWire(40, 100, reverse=True), 25.0)
        o1, o2 = Loft.orientProfiles(p1, p2)
        self.assertLess((o1[0] - o2[0]).Length, (o1[0] - o2[-1]).Length)
        winding = Loft.windingVector(o1, Loft.polylineCenter(o1))
        self.assertGreater(winding.z, 0)

    def test_straight_profiles_align_by_endpoints(self):
        w1 = Part.makePolygon([V(0, 0, 0), V(100, 0, 0)])
        w2 = Part.makePolygon([V(100, 50, 20), V(0, 50, 20)])
        shell = Loft.smLoftShell(w1, w2)
        self.assertEqual(len(shell.Faces), 2)
        self.assertAlmostEqual(shell.Area, 100 * math.hypot(50, 20), places=4)

    def test_all_strategies_build_valid_shell(self):
        w1 = Part.Wire([Part.makeLine(V(-80, -40, 0), V(-50, 0, 0)),
                        _arcWire(50, 0).Edges[0].reversed(),
                        Part.makeLine(V(50, 0, 0), V(80, -40, 0))])
        w2 = _arcWire(30, 120)
        for name in Loft.TRIANGULATION_STRATEGIES:
            with self.subTest(strategy=name):
                shell = Loft.smLoftShell(w1, w2, 25.0, name)
                self.assertEqual(shell.ShapeType, "Shell")
                self.assertTrue(shell.isValid())
                p1 = Loft.discretizeWire(w1, 25.0)
                p2 = Loft.discretizeWire(w2, 25.0)
                self.assertEqual(len(shell.Faces), len(p1) + len(p2) - 2)
                # Shell boundary = two profiles + two end rungs.
                free = [e for e in shell.Edges if len(shell.ancestorsOfType(e, Part.Face)) == 1]
                self.assertEqual(len(free), len(p1) + len(p2))

    def test_default_strategy_is_min_dihedral(self):
        self.assertEqual(Loft.DEFAULT_STRATEGY, "MinDihedral")
        self.assertEqual(list(Loft.TRIANGULATION_STRATEGIES)[:2],
                         ["MinDihedral", "ShortestDiagonal"])

    def test_min_dihedral_first_step_lookahead(self):
        rect = Part.makePolygon([V(50, 0, 0), V(50, 50, 0), V(-50, 50, 0), V(-50, 0, 0)])
        p1, p2 = Loft.orientProfiles(Loft.discretizeWire(rect, 25.0),
                                     Loft.discretizeWire(_arcWire(40, 100), 25.0))
        front = Loft.LoftFront(p1, p2)
        pick = Loft._bestFirstStep(front, front.candidate1(), front.candidate2())
        self.assertIn(pick, (1, 2))
        shell = Loft.makeTriangleShell(Loft.triangulateProfiles(p1, p2))
        self.assertEqual(shell.ShapeType, "Shell")
        self.assertTrue(shell.isValid())

    def test_feature_object(self):
        doc = FreeCAD.newDocument("LoftTest")
        try:
            a = doc.addObject("Part::Feature", "A")
            a.Shape = _arcWire(50, 0)
            b = doc.addObject("Part::Feature", "B")
            b.Shape = _arcWire(30, 100)
            obj = doc.addObject("Part::FeaturePython", "Loft")
            Loft.SMLoft(obj, [a, b])
            doc.recompute()
            self.assertFalse(obj.Shape.isNull())
            self.assertEqual(obj.Shape.ShapeType, "Shell")
            self.assertTrue(obj.Shape.isValid())
            self.assertEqual(obj.MinLineLength.Value, 25.0)
            self.assertEqual(obj.TriangulationMode, "MinDihedral")
            self.assertEqual(obj.StartDirection, V(1, 0, 0))
            a.Shape = Part.makeCircle(50)
            b.Shape = _rectWire(40, 30, 100)
            doc.recompute()
            self.assertEqual(obj.Shape.ShapeType, "Shell")
            self.assertTrue(obj.Shape.isValid())
        finally:
            FreeCAD.closeDocument(doc.Name)

    def test_closed_wire_discretization(self):
        pts = Loft.discretizeWire(Part.Wire(Part.makeCircle(50)), 25.0)
        self.assertEqual(len(pts), int(math.floor(100 * math.pi / 25.0)) + 1)
        self.assertTrue(Loft.isClosedPolyline(pts))
        # A full circle gets at least three chords.
        pts = Loft.discretizeWire(Part.Wire(Part.makeCircle(5)), 25.0)
        self.assertEqual(len(pts), 4)
        rect = Loft.discretizeWire(_rectWire(40, 30, 0), 25.0)
        self.assertEqual(len(rect), 5)
        self.assertTrue(Loft.isClosedPolyline(rect))

    def test_closed_start_vertices(self):
        circle = Part.Wire(Part.makeCircle(50, V(0, 0, 100)))
        rect = _rectWire(40, 30, 0)
        for direction in (V(1, 0, 0), V(0, -1, 0), V(-1, 1, 0)):
            with self.subTest(direction=direction):
                p1, p2 = Loft.prepareProfiles(circle, rect, 25.0, direction)
                center = Loft.polylineCenter(p1)
                best = max((p - center).normalize().dot(V(direction).normalize())
                           for p in p1[:-1])
                self.assertAlmostEqual((p1[0] - center).normalize()
                                       .dot(V(direction).normalize()), best, places=9)
                nearest = min((p - p1[0]).Length for p in p2[:-1])
                self.assertAlmostEqual((p2[0] - p1[0]).Length, nearest, places=9)
                self.assertTrue(Loft.isClosedPolyline(p1))
                self.assertTrue(Loft.isClosedPolyline(p2))

    def test_start_direction_along_axis_is_projected(self):
        circle = Part.Wire(Part.makeCircle(50, V(0, 0, 100)))
        rect = _rectWire(40, 30, 0)
        p1, _ = Loft.prepareProfiles(rect, circle, 25.0, V(0, 0, 1))
        self.assertGreater(p1[0].x, 0)  # falls back to +X

    def test_closed_to_closed_shell(self):
        circle = Part.Wire(Part.makeCircle(50, V(0, 0, 100)))
        rect = _rectWire(40, 30, 0)
        for name in Loft.TRIANGULATION_STRATEGIES:
            with self.subTest(strategy=name):
                shell = Loft.smLoftShell(rect, circle, 25.0, name)
                self.assertEqual(shell.ShapeType, "Shell")
                self.assertTrue(shell.isValid())
                n1 = len(Loft.discretizeWire(rect, 25.0)) - 1
                n2 = len(Loft.discretizeWire(circle, 25.0)) - 1
                self.assertEqual(len(shell.Faces), n1 + n2)
                # Tube: the only free edges are the two profiles.
                free = [e for e in shell.Edges if len(shell.ancestorsOfType(e, Part.Face)) == 1]
                self.assertEqual(len(free), n1 + n2)

    def test_closed_to_open_shell(self):
        circle = Part.Wire(Part.makeCircle(50, V(0, 0, 100)))
        for w1, w2 in ((circle, _arcWire(40, 0)), (_arcWire(40, 0), circle)):
            shell = Loft.smLoftShell(w1, w2)
            self.assertEqual(shell.ShapeType, "Shell")
            self.assertTrue(shell.isValid())
            p1, p2 = Loft.prepareProfiles(w1, w2, 25.0)
            open_, closed = (p2, p1) if Loft.isClosedPolyline(p1) else (p1, p2)
            nearest = min((p - open_[0]).Length for p in closed[:-1])
            self.assertAlmostEqual((closed[0] - open_[0]).Length, nearest, places=9)
