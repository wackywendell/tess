#include "voro++.hh"
#include <boost/python.hpp>
using namespace boost::python;

using namespace voro;

void (container::*container_put)(int, double, double, double) = &container::put;

void (container_poly::*container_poly_put)(int, double, double, double, double) = &container_poly::put;

BOOST_PYTHON_MODULE(pyvoro)
{
    class_<container>("Container", init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_put)
        .def("compute_all_cells", &container::compute_all_cells)
		.def("sum_cell_volumes", &container::sum_cell_volumes)
        ;
    class_<container_poly>("PolyContainer", init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_poly_put)
        .def("compute_all_cells", &container::compute_all_cells)
		.def("sum_cell_volumes", &container::sum_cell_volumes)
        ;
}
