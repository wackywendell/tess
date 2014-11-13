# distutils: language = c++
# distutils: include_dirs = src
# distutils: sources = src/voro++.cc

from __future__ import division

from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cython.operator cimport dereference

cdef extern from "voro++.hh" namespace "voro":
    cdef cppclass container_base:
        pass
    
    cdef cppclass container:
        container(double,double,double,double,double,double,
				int,int,int,cbool,cbool,cbool,int) except +
        cbool compute_cell(voronoicell_neighbor &c,c_loop_all &vl)
        void put(int, double, double, double)
        int total_particles()
    
    cdef cppclass container_poly:
        container_poly(double,double,double,double,double,double,
				int,int,int,cbool,cbool,cbool,int) except +
        cbool compute_cell(voronoicell_neighbor &c, c_loop_all &vl)
        void put(int, double, double, double, double)
        int total_particles()
        
    cdef cppclass voronoicell_neighbor:
        voronoicell()
        double volume()
    
    cdef cppclass c_loop_all:
        c_loop_all(container_base&)
        cbool start()
        cbool inc()
        int pid()
        void pos(double &x, double &y, double &z)
        

cdef class Cell:
    cdef voronoicell_neighbor *thisptr
    def __cinit__(self):
        self.thisptr = new voronoicell_neighbor()
    
    def __dealloc__(self):
        del self.thisptr
    
    def volume(self): return self.thisptr.volume()

cdef class _Container:
    cdef container *thisptr
    def __cinit__(self, double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,cbool xperiodic_,cbool yperiodic_,cbool zperiodic_,int init_mem):
        self.thisptr = new container(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, 
                xperiodic_, yperiodic_, zperiodic_, init_mem)
    
    def __dealloc__(self):
        del self.thisptr
    
    def put(self, int n, double x, double y, double z):
        self.thisptr.put(n,x,y,z)
    
    def get_cells(self):
        cdef container_base *baseptr = (<container_base *>(self.thisptr))
        cdef c_loop_all *vl = new c_loop_all(dereference(baseptr))
        
        cell = Cell()
        
        cdef int vcells_left = self.thisptr.total_particles()
        cdef int id
        cdef double x = 0
        cdef double y = 0
        cdef double z = 0
        
        mylist = [None for _ in range(vcells_left)]
        
        if not vl.start():
            del vl
            raise ValueError("Failed to start loop")
        
        while True:
            if(self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl))):
                cell.id = vl.pid()
                assert(cell.id <= self.thisptr.total_particles())
                mylist[cell.id] = cell
                
                vl.pos(x,y,z)
                cell.pos = (x,y,z)
                
                vcells_left -= 1;
            if not vl.inc(): break
        
        del vl
        
        if vcells_left != 0:
            raise ValueError("Computation failed")

cdef class _ContainerPoly:
    cdef container_poly *thisptr
    def __cinit__(self, double ax_,double bx_,double ay_,double by_,double az_,double bz_,
				int nx_,int ny_,int nz_,cbool xperiodic_,cbool yperiodic_,cbool zperiodic_,int init_mem):
        self.thisptr = new container_poly(ax_, bx_, ay_, by_, az_, bz_, nx_, ny_, nz_, 
                xperiodic_, yperiodic_, zperiodic_, init_mem)
    
    def __dealloc__(self):
        del self.thisptr
    
    def put(self, int n, double x, double y, double z, double r):
        self.thisptr.put(n,x,y,z,r)
    
    def get_cells(self):
        cdef container_base *baseptr = (<container_base *>(self.thisptr))
        cdef c_loop_all *vl = new c_loop_all(dereference(baseptr))
        
        cell = Cell()
        
        cdef int vcells_left = self.thisptr.total_particles()
        cdef int id
        cdef double x = 0
        cdef double y = 0
        cdef double z = 0
        
        mylist = [None for _ in range(vcells_left)]
        
        if not vl.start():
            del vl
            raise ValueError("Failed to start loop")
        
        while True:
            if(self.thisptr.compute_cell(dereference(cell.thisptr), dereference(vl))):
                cell.id = vl.pid()
                assert(cell.id <= self.thisptr.total_particles())
                mylist[cell.id] = cell
                
                vl.pos(x,y,z)
                cell.pos = (x,y,z)
                
                vcells_left -= 1;
            if not vl.inc(): break
        
        del vl
        
        if vcells_left != 0:
            raise ValueError("Computation failed")

class Container(list):
    def __init__(self, points, limits=1.0, periodic=False, radii=None, blocks=None):
        """__init__(self, points, limits=1.0, periodic=False, radii=None, blocks=None)
        
        Get the voronoi cells for a given set of points.
        
        points: an Nx3 iterable of floating point numbers denoting the coordinates.
        limits: the size of the box. May be a number (for a cubic box), or a 3-tuple.
        periodic: a bool or 3-tuple of bools representing wall periodicity
        """
        # make px, py, pz from periodic, whether periodic is a 3-tuple or bool
        try:
            px, py, pz = periodic
        except TypeError:
            px = py = pz = periodic
        px = bool(periodic)
        py = bool(periodic)
        pz = bool(periodic)
        
        # make lx, ly, lz from limits, whether limits is a 3-tuple or float
        try:
            lx, ly, lz = limits
        except TypeError:
            lx = ly = lz = limits
        lx = float(lx)
        ly = float(ly)
        lz = float(lz)
        assert lx > 0
        assert ly > 0
        assert lz > 0
        
        N = len(points)
        
        # make bx, by, bz from blocks, or make it up
        if blocks is None:
            Nthird = pow(N, 1.0/3.0)
            blocks = round(Nthird / lx), round(Nthird / ly), round(Nthird / lz)
        
        try:
            bx, by, bz = blocks
        except TypeError:
            bx = by = bz
        bx = max(int(bx), 1)
        by = max(int(by), 1)
        bz = max(int(bz), 1)
        
        # If we have periodic conditions, we want to get the 'boxed' version of each position.
        # Each coordinate (x,y,z) may or may not be periodic, so we'll deal with them separately.
        def roundedoff(n, l, periodic):
            if periodic:
                return float(n) % l
            else:
                return float(n)
        
        # voro has two types: Container and ContainerPoly. ContainerPoly is for unequal radii,
        # Container is for no-radii.
        # Now we choose the right one.
        if radii is not None:
            assert(len(radii) == len(points))
            self._container = _ContainerPoly(0,lx, 0,ly, 0,lz, # limits
                                bx, by, bz,        # block size
                                px, py, pz,        # periodicity
                                len(points))
            
            for n,(x,y,z),r in zip(range(len(points)), points, radii):
                self._container.put(n, roundedoff(x,lx,px), roundedoff(y,ly,py), roundedoff(z,lz,pz), float(r))
            
            cells = self._container.get_cells()
            list.__init__(self, cells)
            
        else:
            # no radii => use voro._Container
            self._container = _Container(0,lx, 0,ly, 0,lz,         # limits
                                    bx, by, bz,        # block size
                                    px, py, pz,        # periodicity
                                    len(points))
            for n,(x,y,z) in enumerate(points):
                self._container.put(n, roundedoff(x,lx,px), roundedoff(y,ly,py), roundedoff(z,lz,pz))
            
            cells = self._container.get_cells()
            list.__init__(self, cells)
        
        # Sometimes a _Container has calculation issues. That can lead to the following.
        if len(self) != len(points):
            raise ValueError("Points could not be suitably fitted into the given box. \n\
You may want to check that all points are within the box, and none are overlapping. {} / {}".format(len(self), len(points)))
