from setuptools import setup, Extension
import os.path
import sys
import glob

# define the name of the extension to use
extension_name    = 'tess'
extension_version = '0.1'
 
# define the directories to search for include files
# to get this to work, you may need to include the path
# to your boost installation. Mine was in
# '/usr/local/include', hence the corresponding entry.
include_dirs = [ '/usr/local/include', 'src' ]
 
# define the library directories to include any extra
# libraries that may be needed.  The boost::python
# library for me was located in '/usr/local/lib'
library_dirs = [ '/usr/local/lib' ]
 
# define the libraries to link with the boost python library

libraries = [ 'boost_python3' if sys.version_info.major == 3 else 'boost_python' ]
 
# define the source files for the extension
source_files = ['pywrapper.cc', 'cell.cc', 'common.cc', 'container.cc', 'unitcell.cc', 
                'v_compute.cc', 'c_loops.cc', 'v_base.cc', 'wall.cc',
                'pre_container.cc', 'container_prd.cc']
source_files = [os.path.join('src', fname) for fname in source_files]

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# create the extension and add it to the python distribution
setup( name=extension_name, 
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

        ext_modules=[Extension(
            extension_name, source_files, 
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries)])
