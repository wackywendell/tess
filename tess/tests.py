from . import *
from unittest import TestCase

from math import sqrt
import numpy as np

try:
    import scipy
except ImportError:
    scipy = None

class LatticeTest:
    """A basic test for testing a lattice. The basic constants below need to be overwritten,
    as well as the setUp() function."""
    limits = None
    order = None
    
    def setUp(self):
        self.cells = None
    
    def volumes(self, vs):
        found_vs = sorted(set([round(c.volume(), 8) for c in self.cells]))
        self.assertEqual(len(vs), len(found_vs), found_vs)
        
        for v, fv in zip(sorted(vs), found_vs):
            self.assertAlmostEqual(v, fv)
    
    def neighbors(self, ns):
        found_ns = sorted(set([len(c.neighbors()) for c in self.cells]))
        self.assertEqual(len(ns), len(found_ns))
        
        for n, fn in zip(sorted(ns), found_ns):
            self.assertEqual(n, fn)
    
    def test_order(self):
        if scipy is None: return
        if self.order is None: return
        self.assertAlmostEqual(self.order, self.cells.order(l=6, weighted=True, local=True))
        self.assertAlmostEqual(self.order, self.cells.order(l=6, weighted=True, local=False))
        self.assertAlmostEqual(self.order, self.cells.order(l=6, weighted=False, local=True))
        self.assertAlmostEqual(self.order, self.cells.order(l=6, weighted=False, local=False))
    
class CubicLattice(LatticeTest):
    n = 4
    limits = (n,n,n)
    
    def get_points(self, r2 = 0.6):
        nm = self.n
        ptsr = [tuple(np.array((l,m,n))+.5) + ((l+m+n) % 2,)
               for l in range(nm) for m in range(nm) for n in range(nm)]
        ptsr = np.array(ptsr)
        pts = ptsr[:, :3]
        r = [.5 if r < 1 else r2 for r in ptsr[:, 3]]
        return pts, r
        
    def setUp(self):
        pts, _ = self.get_points(r2=0.5)
        self.cells = Container(pts, self.n, periodic=False)
    
    def test_str(self):
        self.assertEqual(str(self.cells[0]), '<Cell 0>')
        self.assertEqual(repr(self.cells[0]), '<Cell 0>')
        self.assertEqual(str(self.cells[1]), '<Cell 1>')
        self.assertEqual(repr(self.cells[1]), '<Cell 1>')
        self.assertEqual(str(self.cells[-1]), '<Cell {0}>'.format(len(self.cells)-1))
        self.assertEqual(repr(self.cells[-1]), '<Cell {0}>'.format(len(self.cells)-1))

class BasicCubic(CubicLattice, TestCase):
    order = 0.35355339059328
    
    def test_volumes(self): self.volumes([1.0])
        
    def test_neighbors(self): self.neighbors([6])
    
    def test_centroid(self):
        cx, cy, cz = self.cells[0].centroid()
        self.assertAlmostEqual(cx, 0.5)
        self.assertAlmostEqual(cy, 0.5)
        self.assertAlmostEqual(cz, 0.5)
    
    def test_face_areas(self):
        areas = self.cells[0].face_areas()
        self.assertEqual(len(areas), 6)
        for a in areas:
            self.assertAlmostEqual(a, 1.0)

class CubicAlternating(CubicLattice, TestCase):
    def setUp(self):
        pts, r = self.get_points(r2=0.6)
        self.cells = Container(pts, self.n, radii=r, periodic=True)
        
    def test_volumes(self): self.volumes([1.295031, 0.704969])

class FCC(LatticeTest, TestCase):
    order = 0.57452425971404164
    n = 4
    limits = np.array((n,n,n))
    
    def points(self):
        x0 = np.array((0,0,0))
        x1 = np.array((0,0.5,0.5))
        x2 = np.array((0.5,0,0.5))
        x3 = np.array((0.5,0.5,0))
        
        dxs = np.array([(i,j,k) for k in range(self.n) for j in range(self.n) for i in range(self.n)])
        xs = np.array([[x+dx for x in (x0, x1, x2, x3)] for dx in dxs])
        return np.concatenate(xs)
    
    def setUp(self):
        self.cells = Container(self.points(), limits=self.limits, periodic=True)
    
    def test_volumes(self): self.volumes([0.25])
    
    def test_neighbors(self): self.neighbors([12])

class HCP(FCC):
    order = 0.48476168522368324
    n = 4
    limits = np.array((n, n*sqrt(3.0)/2, n*sqrt(6.)/3.))
    
    def points(self):
        assert self.n % 2 == 0
        arr = np.array([(
                            2*i + (j%2) + (k%2),
                            sqrt(3.0)*(j + (k%2)/3.),
                            sqrt(24.)*k/3)
                         for k in range(self.n) for j in range(self.n) for i in range(self.n)])
        return arr/2.
    
    def setUp(self):
        self.cells = Container(self.points(), limits=self.limits, periodic=True)
    
    def test_volumes(self): self.volumes([1./sqrt(2.)])
    
    def test_neighbors(self): self.neighbors([12])
