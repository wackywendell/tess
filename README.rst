Tess, a 3D cell-based Voronoi library based on voro++
==============================================

`Original work`_ by Chris H. Rycroft (UC Berkeley / Lawrence Berkeley Laboratory). This repository includes Python bindings, using Boost.Python.
-------------

Voro++ is a software library for carrying out three-dimensional computations
of the Voronoi tessellation. A distinguishing feature of the Voro++ library
is that it carries out cell-based calculations, computing the Voronoi cell
for each particle individually, rather than computing the Voronoi
tessellation as a global network of vertices and edges. It is particularly
well-suited for applications that rely on cell-based statistics, where
features of Voronoi cells (eg. volume, centroid, number of faces) can be
used to analyze a system of particles.

Voro++ comprises of several C++ classes that can be built as a static library
and linked to. A command-line utility is also provided that can analyze text
files of particle configurations and use most of the features of the code.
Numerous examples are provided to demonstrate the library's features and all of
these are discussed in detail on the library website.

.. _Original work: http://math.lbl.gov/voro++/

Tess
======

`This repository`_ is a python wrapper 
around the voro++ library. It is currently not as powerful as the full C++
library, but contains a number of functions.

.. _This repository: https://github.com/wackywendell/tess

Compilation
-----------
To make, run ``python setup.py build_ext`` in either the main directory or the ``src`` directory.

To install, run ``python setup.py install --user`` to install to your home directory (recommended), or
``sudo python setup.py install`` to install to the main setup directory.

Use
---
After installation, use `import tess` to import the module.

Examples of use are `available here`_,
with the notebook being rendered from ``examples/Examples.ipynb``.

.. _available here: http://nbviewer.ipython.org/github/wackywendell/tess/blob/master/examples/Examples.ipynb

Voro++ Copyright Notice
=======================
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


Voro++ Acknowledgments
======================
This work (voro++) was supported by the Director, Office of Science, Computational and
Technology Research, U.S. Department of Energy under Contract No.
DE-AC02-05CH11231.
