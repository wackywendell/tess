from distutils.core import setup, Extension
import os.path
import glob

# define the name of the extension to use
extension_name    = 'pyvoro'
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
libraries = [ 'boost_python' ]
 
# define the source files for the extension
source_files = ['cell.cc', 'common.cc', 'container.cc', 'unitcell.cc', 
                'v_compute.cc', 'c_loops.cc', 'v_base.cc', 'wall.cc',
                'pre_container.cc', 'container_prd.cc']
source_files = [os.path.join('src', fname) for fname in source_files]

print(source_files)

# create the extension and add it to the python distribution
setup( name=extension_name, 
        version=extension_version, 
        ext_modules=[Extension(
            extension_name, source_files, 
            include_dirs=include_dirs,
            library_dirs=library_dirs,
            libraries=libraries)])
