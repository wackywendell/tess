from distutils.core import setup
from distutils.extension import Extension
from distutils.command.sdist import sdist as _sdist
import os.path

try:
    from Cython.Build import cythonize
except ImportError:
    cythonize = None

extension_version = '0.1.2'

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

if cythonize is not None:
    print("Building with Cython.")
    ext = cythonize(Extension("tess._voro",
              sources=["tess/_voro.pyx", "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
              ))
else:
    print("Cython not found, using prebuilt file.")
    ext = [Extension("tess._voro",
              sources=["tess/_voro.cpp", "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
              )]

# Set sdist to make the .cpp file
# from http://stackoverflow.com/a/18418524/4190270
class sdist(_sdist):
    def run(self):
        # this is already imported, but the import might have failed. If so, raise an ImportError now.
        from Cython.Build import cythonize
        
        # Make sure the compiled Cython files in the distribution are up-to-date
        cythonize(['tess/_voro.pyx'])
        _sdist.run(self)

cmdclass = dict(sdist = sdist)

# create the extension and add it to the python distribution
setup( name='tess', 
        version=extension_version, 
        author="Wendell Smith",
        author_email="wackywendell@gmail.com",
        description = ("A module for calculating and analyzing Voronoi tessellations"),
        license = "BSD",
        keywords = "laguerre voronoi tessellation voro++",
        url = "https://github.com/wackywendell/tess",
        long_description=read('README.rst'),
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Topic :: Scientific/Engineering :: Mathematics",
            "Topic :: Scientific/Engineering :: Physics",
            "Intended Audience :: Science/Research",
            "Operating System :: POSIX :: Linux",
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: BSD License",
        ],
        
        packages=['tess'],
        package_dir={"tess": "tess"},
        
        ext_modules = ext,
        cmdclass=cmdclass,
        )
