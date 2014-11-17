Quick Start
-----------

Installation
~~~~~~~~~~~~

To install, use ``pip`` (or ``easy_install``)::

    pip install --user tess

Or to install from Github_::
    
    pip install --user git+git://github.com/wackywendell/tess@master
    
.. _github: https://www.github.com/wackywendell/tess

Usage
~~~~~

The first step is to create a :class:`tess.Cell`::

    >>> cntr = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)

A container is a list of :class:`tess.Cell` objects, representing Voronoi ells:
    
    >>> [round(v.volume(), 3) for v in cntr]
    [13.5, 13.5]

:class:`tess.Cell` objects have many methods. Here are a few::

    >>> [v.pos for v in cntr]
    [(1.0, 1.0, 1.0), (2.0, 2.0, 2.0)]
    
    >>> [v.centroid() for v in cntr]
    [(1.09375, 1.09375, 1.09375), (1.90625, 1.90625, 1.90625)]
    
    >>> [v.neighbors() for v in cntr]
    [[-5, -2, -3, -1, -4, 1, -6], [0, -3, -6, -4, -5, -2, -1]]
    
    >>> [v.face_areas() for v in cntr]
    [[7.875, 1.125, 7.875, 7.875, 1.125, 11.691342951089922, 1.125],
     [11.691342951089922, 1.125, 7.875, 7.875, 1.125, 7.875, 1.125]]
    
    >>> [v.normals() for v in cntr]
    [[(0.0, 0.0, -1.0),
      (1.0, 0.0, 0.0),
      (0.0, -1.0, 0.0),
      (-1.0, 0.0, 0.0),
      (0.0, 1.0, 0.0),
      (0.5773502691896257, 0.5773502691896257, 0.5773502691896257),
      (0.0, 0.0, 1.0)],
     [(-0.5773502691896257, -0.5773502691896257, -0.5773502691896257),
      (-0.0, -1.0, -0.0),
      (0.0, 0.0, 1.0),
      (0.0, 1.0, -0.0),
      (0.0, 0.0, -1.0),
      (1.0, 0.0, -0.0),
      (-1.0, -0.0, -0.0)]]
      
See the Reference_ for more methods, or just use a Python interpreter or IPython notebook to find
them on your own!

.. _Reference: api.html
