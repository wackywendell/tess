"""
Tess
****

A 3D cell-based Voronoi library based on voro++
-----------------------------------------------

`Code available`_ on Github.

`Full documentation available`_ at Read the Docs.

.. _Code available: https://github.com/wackywendell/tess

.. _Full documentation available: https://tess.readthedocs.org

Description
-----------

Tess is a library to calculate Voronoi (and Laguerre) tessellations in 3D and analyze their
structure. The tessellation is calculated as a `list` of `Cell` objects, each of which
can give information on its volume, centroid, number of faces, surface area, etc. The library is
made with packings of spherical particles in mind, possibly with variable sizes.

"""
import importlib
from setuptools import setup, Extension
from setuptools.command.sdist import sdist


# setuptools transparently fallback to using the cpp file if cython is not available. See:
# https://setuptools.pypa.io/en/latest/userguide/distribution.html#distributing-extensions-compiled-with-cython
extension = Extension(
            "tess._voro",
            sources=["tess/_voro.pyx", "src/voro++.cc"],
            include_dirs=["src"],
            language="c++",
)


class cython_sdist(sdist):
    # Set sdist to make the .cpp file
    # from http://stackoverflow.com/a/18418524/4190270

    def run(self):
        # this is already imported, but the import might have failed.
        # If so, raise an ImportError now.
        from Cython.Build import cythonize
        # Make sure the compiled Cython files in the distribution are up-to-date
        cythonize(["tess/_voro.pyx"])
        sdist.run(self)


# create the extension and add it to the python distribution
setup(
    name="tess",
    author="Wendell Smith",
    author_email="wackywendell@gmail.com",
    description=("A module for calculating and analyzing Voronoi tessellations"),
    license="BSD",
    keywords="laguerre voronoi tessellation voro++",
    url="https://tess.readthedocs.org",
    long_description=__doc__.strip(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Scientific/Engineering :: Mathematics",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
    ],
    packages=["tess"],
    package_dir={"tess": "tess"},
    cmdclass={'sdist': cython_sdist},
    ext_modules=[extension],
    use_scm_version=True,
    setup_requires=["setuptools>=40.9.0", "wheel", "setuptools_scm"],
    extras_require={
        'tests': ["numpy", "scipy", "sphinx", "sphinx_rtd_theme"],
    },
)
