"""
This is a library to calculate Voronoi cells and access their information.

Basic Process
~~~~~~~~~~~~~

  - Create a :class:`Container` object, using information about your system. 
      - a  :class:`Container` is a `list` of :class:`Cell` objects
  - Access the :class:`Cell` methods to get information about them

Example
~~~~~~~

    >>> from tess import Container
    >>> c = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)
    >>> [round(v.volume(), 3) for v in c]
    [13.5, 13.5]
"""

from ._voro import Container as _Container, ContainerPoly as _ContainerPoly, Cell


class Container(list):
    r"""A container (`list`) of Voronoi cells.
    
    This is the main entry point into the :mod:`tess` module. After creation, this will be a `list` 
    of :class:`Cell` objects.
    
    The :class:`Container` must be rectilinear, and can have solid boundary conditions, periodic 
    boundary conditions, or a mix of the two.
    
    >>> from tess import Container
    >>> c = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)
    >>> [round(v.volume(), 3) for v in c]
    [13.5, 13.5]
    
    Parameters
    ----------
    points : iterable of iterable of `float`
        The coordinates of the points, size Nx3.
    limits : `float`, 3-tuple of float, or two 3-tuples of float
        The box limits. 
        If given a float `L`, then the box limits are [0, 0, 0] to [L, L, L].
        If given a 3-tuple `(Lx, Ly, Lz)`, limits are [0, 0, 0] to [Lx, Ly, Lz].
        If given two 3-tuples `(x0, y0, z0), (x1, y1, z1)`, limits are [x0, y0, z0] to [x1, y1, z1].
    periodic : `bool` or 3-tuple of `bool`, optional
        Periodicity of the x, y, and z walls
    radii : iterable of `float`, optional
        for unequally sized particles, for generating a Laguerre transformation.
    
    Returns
    -------
    Container
        A `list` of :class:`Cell` objects
    
    Notes
    -----
    *Voronoi Tesselation*
    
    A point :math:`\vec x` is part of a Voronoi cell :math:`i` with nucleus :math:`\vec{r}_i` iff
    
    .. math:: 
        \left|\vec{x}-\vec{r}_{i}\right|^{2}<\left|\vec{x}-\vec{r}_{j}\right|^{2} \forall j\neq i
    
    
    *Laguerre Tesselation*, also known as *Radical Voronoi Tesselation*
    
    A point :math:`\vec x` is part of a Laguerre cell :math:`i` with nucleus :math:`\vec{r}_i` and
    radius :math:`R_i` iff
    
    .. math::
        \left|\vec{x}-\vec{r}_{i}\right|^{2}-R_{i}^{2}<\left|\vec{x}-\vec{r}_{j}\right|^{2} -
        R_{j}^{2}\forall j\neq i
    """
    
    def __init__(self, points, limits=1.0, periodic=False, radii=None, blocks=None):
        """Get the voronoi cells for a given set of points."""
        # make px, py, pz from periodic, whether periodic is a 3-tuple or bool
        try:
            px, py, pz = periodic
        except TypeError:
            px = py = pz = periodic
        px = bool(periodic)
        py = bool(periodic)
        pz = bool(periodic)
        
        # make lx, ly, lz from limits, whether limits is a 3-tuple or float
        # lx0, ly0, lz0 are the lower limits, lx, ly, lz are the upper limits
        # If `periodic`, then only the relative difference matters
        try:
            l0, far_limits = limits
            lx0, ly0, lz0 = l0
            lx, ly, lz = far_limits
        except ValueError:
            lx0 = ly0 = lz0 = 0.0
            lx, ly, lz = limits
        except TypeError:
            lx0 = ly0 = lz0 = 0.0
            lx = ly = lz = limits
        lx0 = float(lx0)
        ly0 = float(ly0)
        ly0 = float(lz0)
        lx = float(lx)
        ly = float(ly)
        lz = float(lz)
        assert lx > lx0
        assert ly > ly0
        assert lz > lz0
        
        # The actual lengths on each side
        Lx, Ly, Lz = lx - lx0, ly - ly0, lz - lz0
        
        N = len(points)
        
        # make bx, by, bz from blocks, or make it up
        if blocks is None:
            V = Lx * Ly * Lz
            Nthird = pow(N/V, 1.0/3.0)
            blocks = round(Nthird * Lx), round(Nthird * Ly), round(Nthird * Lz)
        
        try:
            bx, by, bz = blocks
        except TypeError:
            bx = by = bz = blocks
        bx = max(int(bx), 1)
        by = max(int(by), 1)
        bz = max(int(bz), 1)
        self.blocks = bx, by, bz
        
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
            self._container = _ContainerPoly(lx0,lx, ly0,ly, lz0,lz, # limits
                                bx, by, bz,        # block size
                                px, py, pz,        # periodicity
                                len(points))
            
            for n,(x,y,z),r in zip(range(len(points)), points, radii):
                self._container.put(n, roundedoff(x,lx,px), roundedoff(y,ly,py), roundedoff(z,lz,pz), float(r))
        else:
            # no radii => use voro._Container
            self._container = _Container(lx0,lx, ly0,ly, lz0,lz,         # limits
                                    bx, by, bz,        # block size
                                    px, py, pz,        # periodicity
                                    len(points))
            for n,(x,y,z) in enumerate(points):
                self._container.put(n, roundedoff(x,Lx,px), roundedoff(y,Ly,py), roundedoff(z,Lz,pz))
        
        cells = self._container.get_cells()
        list.__init__(self, cells)
        
        # Sometimes a _Container has calculation issues. That can lead to the following.
        if len(self) != len(points):
            raise ValueError("Points could not be suitably fitted into the given box. \n\
You may want to check that all points are within the box, and none are overlapping. {} / {}".format(len(self), len(points)))
    
    def get_widths(self):
        "Get the size of the box."
        return self._container.get_widths()
        
    def _get_bond_normals(self):
        """Returns a generator of [(dx,dy,dz,A) for each bond] for each cell.
        
        (dx,dy,dz) is the normal, and A is the area of the voronoi face."""
        return [
            [(x,y,z,A) for (x,y,z), A in zip(vc.normals(), vc.face_areas())]
            for vc in self
        ]

    def order(self, l=6, local=False, weighted=True):
        r"""Returns crystalline order parameter :math:`Q_l` (such as :math:`Q_6`).
        
        Requires numpy and scipy.
        
        Parameters
        ----------
        l : int, optional
            Defines which :math:`Q_l` you want (6 is standard, for detecting hexagonal lattices)
        local : bool, optional
            Calculate Local :math:`Q_6` (true) or Global :math:`Q_6`
        weighted : bool, optional
            Whether or not to weight by area the faces of each polygonal side
        
        Notes
        -----
        
        For ``local=False``, this calculates
        
        .. math:: 
            Q_l = \sqrt{\frac{4 \pi}{2 l + 1}\sum_{m=-l}^{l} \left| 
                    \sum_{i=1}^{N_b} w_i Y_{lm}\left(\theta_i, \phi_i \right) \right|^2}
        
        where:
        
        :math:`N_b` is the number of bonds
        
        :math:`\theta_i` and :math:`\phi_i` are the angles of each bond :math:`i`, in spherical 
        coordinates
        
        :math:`Y_{lm}\left(\theta_i, \phi_i \right)` is the spherical harmonic function
        
        :math:`w_i` is the weighting factor, either proportional to the area (for `weighted`)
        or all equal (:math:`\frac{1}{N_b}`)
        
        For ``local=True``, this calculates
        
        .. math:: 
            Q_{l,\mathrm{local}} = \sum_{j=1}^N \sqrt{\frac{4 \pi}{2 l + 1}\sum_{m=-l}^{l} \left| 
                    \sum_{i=1}^{n_b^j} w_i Y_{lm}\left(\theta_i, \phi_i \right) \right|^2}
                    
        where variables are as above, and each *cell* is weighted equally but each *bond* for each 
        cell is weighted: :math:`\sum_{i=1}^{n_b^j} w_i = 1`
        
        Returns
        -------
        
        float
        
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
    """Converts 3D cartesian coordinates to the angular portion of spherical coordinates, (theta, phi).
    
    Requires numpy.
    
    Parameters
    ----------
    
    xyz: array-like, Nx3
        
        Column 0: the "elevation" angle, :math:`0` to :math:`\pi`
        
        Column 1: the "azimuthal" angle, :math:`0` to :math:`2\pi`
    
    Returns
    -------
    array, Nx2
    """
    import numpy as np
    ptsnew = np.zeros((xyz.shape[0],2))
    xy = xyz[:,0]**2 + xyz[:,1]**2
    ptsnew[:,0] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
    ptsnew[:,1] = np.arctan2(xyz[:,1], xyz[:,0])
    return ptsnew


def orderQ(l, xyz, weights=1):
    r"""Returns :math:`Q_l`, for a given l (int) and a set of Cartesian coordinates xyz.
    
    Requires numpy and scipy.
    
    For global :math:`Q_6`, use :math:`l=6`, and pass xyz of all the bonds.
    
    For local :math:`Q_6`, use :math:`l=6`, and the bonds have to be averaged slightly differently.
    
    Parameters
    ----------
    l : int
        The order of :math:`Q_l`
    xyz : array-like Nx3
        The bond vectors :math:`\vec r_j - \vec r_i`
    weights : array-like, optional
        How to weight the bonds; weighting by Voronoi face area is common.
    
    Notes
    -----
    This calculates
    
    .. math:: 
        Q_l = \sqrt{\frac{4 \pi}{2 l + 1}\sum_{m=-l}^{l} \left| 
                \sum_{i=1}^{N_b} w_i Y_{lm}\left(\theta_i, \phi_i \right) \right|^2}
    
    where:
    
    :math:`N_b` is the number of bonds
    
    :math:`\theta_i` and :math:`\phi_i` are the angles of each bond :math:`i`, in spherical 
    coordinates
    
    :math:`Y_{lm}\left(\theta_i, \phi_i \right)` is the spherical harmonic function
    
    :math:`w_i` are the `weights`, defaulting to uniform: (:math:`\frac{1}{N_b}`)
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
