from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import os.path

extension_version = '0.1'

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
    
ext = Extension("tess._voro",
              sources=["tess/_voro.pyx", "src/voro++.cc"],
              include_dirs=["src"],
              language="c++",
              )
extensions = cythonize(ext)

# create the extension and add it to the python distribution
setup( name='tess', 
        version=extension_version, 
        author="Wendell Smith",
        author_email="wackywendell@gmail.com",
        description = ("A odule for calculating and analyzing Voronoi tesselations"),
        license = "BSD",
        keywords = "laguerre voronoi tesselation voro++",
        url = "https://github.com/wackywendell/voronoia",
        long_description=read('README.md'),
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
        
        ext_modules = extensions)
