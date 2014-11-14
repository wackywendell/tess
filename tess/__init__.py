from ._voro import Container as _Container, ContainerPoly as _ContainerPoly, Cell

class Container(list):
    def __init__(self, points, limits=1.0, periodic=False, radii=None, blocks=None):
        """Get the voronoi cells for a given set of points.
        
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
    
    def get_widths(self):
        return self._container.get_widths()
        
    def _get_bond_normals(self):
        """Returns a generator of [(dx,dy,dz,A) for each bond] for each cell.
        
        (dx,dy,dz) is the normal, and A is the area of the voronoi face."""
        return [
            [(x,y,z,A) for (x,y,z), A in zip(vc.normals(), vc.face_areas())]
            for vc in self
        ]

    def order(self, l=6, local=False, weighted=True):
        """Returns crystalline order parameter.
        
        Requires numpy and scipy.
        
        l is the which Q you want (6 is standard, for detecting hexagonal lattices)
        local is whether to calculate the Local Q6 or Global Q6
        weighted is whether or not to weight by area of the faces of each polygonal side
        """
        import numpy as np
        
        bonds = self._get_bond_normals()
        Nb = np.sum([len(blst) for blst in bonds])
        if not local:
            bonds = [np.concatenate(bonds)]
        
        Qs = []
        for blst in bonds:
            blst = np.array(blst)
            xyz, A = blst[:,:3], blst[:,3]
            Q = orderQ(l, xyz, weights=(A if weighted else 1))
            if local:
                Q = Q * float(np.shape(xyz)[0])
            Qs.append(Q)
        if local:
            return np.sum(Qs) / Nb
        else:
            assert(len(Qs) == 1)
        return np.mean(Qs)

def cart_to_spher(xyz):
    """Takes Nx3 matrix of Cartesian coordinates, and converts them to (theta, phi).
    
    Requires numpy.
    
    Returns an Nx2 matrix.
    
    Column 0: the "elevation" angle, 0 to pi
    Column 1: the "azimuthal angle, 0 to 2*pi
    """
    import numpy as np
    ptsnew = np.zeros((xyz.shape[0],2))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    ptsnew[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew

def orderQ(l, xyz, weights=1):
    """Returns Ql, for a given l (int) and a set of Cartesian coordinates xyz.
    
    Requires numpy and scipy.
    
    For global Q6, use l=6, and pass xyz of all the bonds.
    For local Q6, use l=6, and pass xyz of the bonds of each atom separately, and then average the results.
    
    To weight by Voronoi neighbors, pass weights=(face areas).
    """
    import numpy as np
    from scipy.special import sph_harm
    
    theta, phi = cart_to_spher(xyz).T
    weights = np.ones(xyz.shape[0]) * weights
    weights /= np.sum(weights)
    mmeans = np.zeros(2*l+1, dtype=float)
    for m in range(-l, l+1):
        sph_weighted = sph_harm(m, l, phi, theta).dot(weights) # Average of Y_{6m}
        mmeans[m] = abs(sph_weighted)**2
    return np.sqrt(4*np.pi/(2*l+1) * np.sum(mmeans))
