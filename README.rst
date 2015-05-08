Tess
****

A 3D cell-based Voronoi library based on voro++
-----------------------------------------------

.. image:: https://travis-ci.org/wackywendell/tess.svg?branch=master
    :target: https://travis-ci.org/wackywendell/tess
    :alt: Build Status


.. image:: https://readthedocs.org/projects/tess/badge/?version=latest
    :target: https://readthedocs.org/projects/tess/?badge=latest
    :alt: Documentation Status

This library includes Python bindings, using Cython.

`Code available`_ on Github.

`Documentation available`_ at Read the Docs. 

.. _Code available: https://github.com/wackywendell/tess

.. _Documentation available: https://tess.readthedocs.org

Description
-----------

Tess is a library to calculate Voronoi (and Laguerre) tessellations in 3D and analyze their
structure. The tessellation is calculated as a list of :class:`~tess.Cell` objects, each of which
can give information on its volume, centroid, number of faces, surface area, etc. The library is 
made with packings of spherical particles in mind, possibly with variable sizes. 

voro++
~~~~~~

The Tess library is a set of Python bindings to the Voro++ library. Voro++ provides all the 
algorithms, and Tess provides an easy to use interface to the voro++ library for Python, using
Cython to do so. 

`Original work`_ on voro++ by Chris H. Rycroft (UC Berkeley / Lawrence Berkeley Laboratory). 

.. _Original work: http://math.lbl.gov/voro++/



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

The first step is to create a :class:`~tess.Container`::
    
    >>> from tess import Container
    >>> cntr = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)

A container is a list of :class:`~tess.Cell` objects, representing Voronoi cells:
    
    .. testsetup:: *
        
        from tess import Container
        cntr = Container([[1,1,1], [2,2,2]], limits=(3,3,3), periodic=False)
    
    >>> [round(v.volume(), 3) for v in cntr]
    [13.5, 13.5]

:class:`~tess.Cell` objects have many methods. Here are a few::

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


Voro++ Copyright And Acknowledgments
------------------------------------

Copyright Notice
~~~~~~~~~~~~~~~~

Voro++ Copyright (c) 2008, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required
approvals from the U.S. Dept. of Energy). All rights reserved.

If you have questions about your rights to use or distribute this software,
please contact Berkeley Lab's Technology Transfer Department at TTD@lbl.gov.

NOTICE. This software was developed under partial funding from the U.S.
Department of Energy. As such, the U.S. Government has been granted for itself
and others acting on its behalf a paid-up, nonexclusive, irrevocable, worldwide
license in the Software to reproduce, prepare derivative works, and perform
publicly and display publicly. Beginning five (5) years after the date
permission to assert copyright is obtained from the U.S. Department of Energy,
and subject to any subsequent five (5) year renewals, the U.S. Government is
granted for itself and others acting on its behalf a paid-up, nonexclusive,
irrevocable, worldwide license in the Software to reproduce, prepare derivative
works, distribute copies to the public, perform publicly and display publicly,
and to permit others to do so.


Acknowledgments
~~~~~~~~~~~~~~~
This work (voro++) was supported by the Director, Office of Science, Computational and
Technology Research, U.S. Department of Energy under Contract No.
DE-AC02-05CH11231.
