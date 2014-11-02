#include "voro++.hh"
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
namespace py = boost::python;

using namespace voro;

#define sptr boost::shared_ptr

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

py::tuple get_first_cell(container& cntr){
    c_loop_all vl(cntr);
    voronoicell c;
    
    if(vl.start()) do if(cntr.compute_cell(c,vl)) {
        return py::make_tuple(vl.pid(), c.volume(), c);
    } while(vl.inc());
    
    c.init_octahedron(1);
    return py::make_tuple(-1, 0, c);
}

voronoicell get_just_cell(container& cntr){
    c_loop_all vl(cntr);
    voronoicell c;
    
    if(vl.start()) do if(cntr.compute_cell(c,vl)) {
        return c;
    } while(vl.inc());
    
    c.init_octahedron(1);
    return c;
}

voronoicell get_one_cell(container& cntr){
    voronoicell c;
    c.init_octahedron(1);
    return c;
}

voronoicell* octahedron(){
    voronoicell* c = new voronoicell();
    c->init_octahedron(1);
    return c;
}

sptr<voronoicell> octahedron_s(){
    sptr<voronoicell> c = sptr<voronoicell>(new voronoicell);
    c->init_octahedron(1);
    return c;
}

py::list get_cells(container& cntr){
    c_loop_all vl(cntr);
    sptr<voronoicell> c = sptr<voronoicell>(new voronoicell());
    py::list my_list = py::list();
    
    if(vl.start()) do if(cntr.compute_cell(*c,vl)) {
        my_list.append(c);
        c = sptr<voronoicell>(new voronoicell());
    } while(vl.inc());
            
    return my_list;
}

py::list get_volumes(container& cntr){
    c_loop_all vl(cntr);
    voronoicell c;
    py::list my_list = py::list();
    
    if(vl.start()) do if(cntr.compute_cell(c,vl)) {
        my_list.append(c.volume());
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
// Instantiate Templates
//~ 
//~ namespace boost {
  //~ namespace python {
    //~ template <class T> struct shared_ptr_voronoicell< shared_ptr<voronoicell> >
    //~ { typedef T type; };
  //~ }}

////////////////////////////////////////////////////////////////////////////////////////////////////
// Make the Python Module

BOOST_PYTHON_MODULE(pyvoro)
{
    py::def("octahedron", octahedron, py::return_value_policy<py::manage_new_object>());
    py::def("octahedrons", octahedron_s);
    
    py::class_<container>("Container", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_put)
        .def("add", container_add)
        .def("compute_all_cells", &container::compute_all_cells)
        .def("get_cells", get_cells)
        .def("get_first_cell", get_first_cell)
        .def("get_just_cell", get_just_cell)
        .def("get_one_cell", get_one_cell)
		.def("sum_cell_volumes", &container::sum_cell_volumes)
		.def("cell_volumes", get_volumes)
        ;
    py::class_<container_poly>("PolyContainer", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_poly_put)
        .def("compute_all_cells", &container_poly::compute_all_cells)
		.def("sum_cell_volumes", &container_poly::sum_cell_volumes)
        ;
        
    py::class_<voronoicell_base>("VoronoiCellBase", py::init<>())
        .def("volume", &voronoicell_base::volume)
        .def("max_radius_squared", &voronoicell_base::max_radius_squared)
        .def("total_edge_distance", &voronoicell_base::total_edge_distance)
        .def("surface_area", &voronoicell_base::surface_area)
        .def("centroid", cell_centroid)
        .def("number_of_faces", &voronoicell_base::number_of_faces)
    ;
    
    py::class_<voronoicell, sptr<voronoicell>, py::bases<voronoicell_base> >("VoronoiCell", py::init<>())
        .def("init", &voronoicell::init)
        .def("init_tetrahedron", &voronoicell::init_tetrahedron)
        .def("init_octahedron", &voronoicell::init_octahedron)
    ;
    
}
