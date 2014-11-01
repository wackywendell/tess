#include "voro++.hh"
#include <boost/python.hpp>
namespace py = boost::python;

using namespace voro;

// Helper function for converting std::vector to a python list
// Taken from http://stackoverflow.com/questions/6157409/stdvector-to-boostpythonlist
template<class T>
py::list std_vector_to_py_list(const std::vector<T>& v)
{
    py::object get_iter = py::iterator<std::vector<T> >();
    py::object iter = get_iter(v);
    py::list l(iter);
    return l;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for resolving Python calls

void (container::*container_put)(int, double, double, double) = &container::put;

void (container_poly::*container_poly_put)(int, double, double, double, double) = &container_poly::put;

py::list get_cells(container& cntr){
    c_loop_all vl(cntr);
    voronoicell c;
    py::list my_list = py::list();
    
    if(vl.start()) do if(cntr.compute_cell(c,vl)) {
        my_list.append(c);
    } while(vl.inc());
            
    return my_list;
}

py::tuple cell_centroid(voronoicell_base& vc){
    double x,y,z;
    vc.centroid(x,y,z);
    return py::make_tuple(x,y,z);
}

// Adapted from http://stackoverflow.com/questions/3761391/boostpython-python-list-to-stdvector
void container_add(container& cntr, py::list& lst){
    for (int i = 0; i < py::len(lst); ++i){
        py::tuple curpos = py::extract<py::tuple>(lst[i]);
        double x = py::extract<double>(curpos[0]);
        double y = py::extract<double>(curpos[1]);
        double z = py::extract<double>(curpos[2]);
        cntr.put(i, x, y, z);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Make the Python Module

BOOST_PYTHON_MODULE(pyvoro)
{
    py::class_<container>("Container", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_put)
        .def("add", container_add)
        .def("compute_all_cells", &container::compute_all_cells)
        .def("get_cells", get_cells)
		.def("sum_cell_volumes", &container::sum_cell_volumes)
        ;
    py::class_<container_poly>("PolyContainer", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_poly_put)
        .def("compute_all_cells", &container_poly::compute_all_cells)
		.def("sum_cell_volumes", &container_poly::sum_cell_volumes)
        ;
        
    py::class_<voronoicell_base>("VoronoiCellBase", py::no_init)
        .def("volume", &voronoicell_base::volume)
        .def("max_radius_squared", &voronoicell_base::max_radius_squared)
        .def("total_edge_distance", &voronoicell_base::total_edge_distance)
        .def("surface_area", &voronoicell_base::surface_area)
        .def("centroid", cell_centroid)
        .def("number_of_faces", &voronoicell_base::number_of_faces)
    ;
    
    py::class_<voronoicell, py::bases<voronoicell_base> >("VoronoiCell", py::no_init)
    ;
    
}
