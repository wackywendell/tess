# Voro++ makefile
#
# Author : Chris H. Rycroft (LBL / UC Berkeley)
# Email  : chr@alum.mit.edu
# Date   : August 30th 2011

# Load the common configuration file
include ../config.mk

# List of the common source files
objs=cell.o common.o container.o unitcell.o v_compute.o c_loops.o \
     v_base.o wall.o pre_container.o container_prd.o
src=$(patsubst %.o,%.cc,$(objs))

# Rules for the python module
PYTHON_INCLUDE = `python-config --includes`
PYTHON_LIB = `python-config --libs`
PYTHON_LD = `python-config --ldflags`

# location of the Boost Python include files and library
BOOST_INC = /usr/include
BOOST_LIB = /usr/lib

# Makefile rules
all: libvoro++.a voro++

depend:
	$(CXX) -MM $(src) >Makefile.dep

include Makefile.dep

libvoro++.a: $(objs)
	rm -f libvoro++.a
	ar rs libvoro++.a $^

voro++: libvoro++.a cmd_line.cc
	$(CXX) $(CFLAGS) -L. -o voro++ cmd_line.cc -lvoro++

%.o: %.cc Makefile
	$(CXX) $(CFLAGS) -fPIC -c $<
	
pyvoro.o: pyvoro.cc libvoro++.a
	g++ $(PYTHON_INCLUDE) -I$(BOOST_INC) -fPIC -L. -lvoro++ -c pyvoro.cc
	
pyvoro.so: pyvoro.o
	g++ -shared -Wl,--export-dynamic pyvoro.o -L. -lvoro++ -L$(BOOST_LIB) -lboost_python3 $(PYTHON_LIB) $(PYTHON_LD) -o pyvoro.so

help: Doxyfile $(SOURCE)
	doxygen Doxyfile

clean:
	rm -f $(objs) voro++ libvoro++.a pyvoro.o pyvoro.so

.PHONY: all help execs depend
