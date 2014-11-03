#include "voro++.hh"
#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
namespace py = boost::python;

using namespace voro;

#define sptr boost::shared_ptr

// Helper function for converting std::vector to a python list
// Taken from http://stackoverflow.com/questions/6157409/stdvector-to-boostpythonlist
template<class T>
py::list std_vector_to_py_list(const std::vector<T>& v){
    py::list l;
    typename std::vector<T>::const_iterator it;
    for (it = v.begin(); it != v.end(); ++it)
    l.append(*it);   
    return l;  
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Structs for storing and accessing cell cata

// For holding info about a neighbor.
// To be used with c_loop: you should set the id and x,y,z at the same time as you compute the cell:
// 
//     c_loop_all vl(cntr);
//     sptr<NeighborCell> cell = sptr<NeighborCell>(new NeighborCell);
//     if(vl.start()) do if(cntr.compute_cell(cell->cell,vl)) {
//         cell->init(vl);
//         ...
//    } while(vl.inc());     

class NeighborCell {
    public:
        int id;
        double x,y,z;
        voronoicell_neighbor cell;
        
    public:
        NeighborCell() : id(-1), x(0), y(0), z(0){};
        inline void init(c_loop_base &cl){id = cl.pid(); cl.pos(x,y,z);}
        
        py::tuple pos(){return py::make_tuple(x,y,z);}
        py::tuple centroid(){
            double cx,cy,cz;
            cell.centroid(cx,cy,cz);
            return py::make_tuple(cx+x,cy+y,cz+z);
        }
    
        inline double volume(){return cell.volume();}
        inline double max_radius_squared(){return cell.max_radius_squared();}
        inline double total_edge_distance(){return cell.total_edge_distance();}
        inline double surface_area(){return cell.surface_area();}
        inline double number_of_faces(){return cell.number_of_faces();}
        
        py::list neighbors(){
            std::vector<int> ns;
            cell.neighbors(ns);
            return std_vector_to_py_list(ns);
        }
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for resolving Python calls

void (container::*container_put)(int, double, double, double) = &container::put;

void (container_poly::*container_poly_put)(int, double, double, double, double) = &container_poly::put;

template <class T>
py::list get_cells(T& cntr){
    c_loop_all vl(cntr);
    sptr<NeighborCell> c = sptr<NeighborCell>(new NeighborCell());
    py::list my_list = py::list();
    
    if(vl.start()) do if(cntr.compute_cell(c->cell,vl)) {
        c->init(vl);
        my_list.append(c);
        c = sptr<NeighborCell>(new NeighborCell());
    } while(vl.inc());
            
    return my_list;
}

inline py::list get_cells_container(container& cntr){return get_cells(cntr);}
inline py::list get_cells_container_poly(container_poly& cntr){return get_cells(cntr);}

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

void container_poly_add(container_poly& cntr, py::list& lst){
    for (int i = 0; i < py::len(lst); ++i){
        py::tuple curpos = py::extract<py::tuple>(lst[i]);
        double x = py::extract<double>(curpos[0]);
        double y = py::extract<double>(curpos[1]);
        double z = py::extract<double>(curpos[2]);
        double r = py::extract<double>(curpos[3]);
        cntr.put(i, x, y, z, r);
    }
}

py::list vc_neighbors(voronoicell_neighbor &vc){
    std::vector<int> ns;
    vc.neighbors(ns);
    return std_vector_to_py_list(ns);
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
    py::class_<container>("Container", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_put)
        .def("add", container_add)
        .def("get_cells", get_cells_container)
 		.def("sum_cell_volumes", &container::sum_cell_volumes)
		.def("cell_volumes", get_volumes)
        ;
    py::class_<container_poly>("ContainerPoly", py::init<double,double,double,double,double,double,
				int,int,int,bool,bool,bool,int>())
        .def("put", container_poly_put)
        .def("add", container_poly_add)
        .def("get_cells", get_cells_container_poly)
        .def("compute_all_cells", &container_poly::compute_all_cells)
		.def("sum_cell_volumes", &container_poly::sum_cell_volumes)
        ;
        
        py::class_<NeighborCell, boost::noncopyable, sptr<NeighborCell> >("Cell", py::no_init)
            .def_readonly("x", &NeighborCell::x)
            .def_readonly("y", &NeighborCell::y)
            .def_readonly("z", &NeighborCell::z)
            .def_readonly("id", &NeighborCell::id)
            .def("pos", &NeighborCell::pos)
            .def("volume", &NeighborCell::volume)
            .def("neighbors", &NeighborCell::neighbors)
            .def("max_radius_squared", &NeighborCell::max_radius_squared)
            .def("total_edge_distance", &NeighborCell::total_edge_distance)
            .def("surface_area", &NeighborCell::surface_area)
            .def("centroid", &NeighborCell::centroid)
            .def("number_of_faces", &NeighborCell::number_of_faces)
        ;
        
    //~ py::class_<voronoicell_base>("VoronoiCellBase", py::init<>())
        //~ .def("volume", &voronoicell_base::volume)
        //~ .def("max_radius_squared", &voronoicell_base::max_radius_squared)
        //~ .def("total_edge_distance", &voronoicell_base::total_edge_distance)
        //~ .def("surface_area", &voronoicell_base::surface_area)
        //~ .def("centroid", cell_centroid)
        //~ .def("number_of_faces", &voronoicell_base::number_of_faces)
    //~ ;
    //~ 
    //~ py::class_<voronoicell, boost::noncopyable, py::bases<voronoicell_base>, sptr<voronoicell> >(
                //~ "VoronoiCell", py::init<>())
        //~ .def("init", &voronoicell::init)
        //~ .def("init_tetrahedron", &voronoicell::init_tetrahedron)
        //~ .def("init_octahedron", &voronoicell::init_octahedron)
    //~ ;
    //~ 
    //~ py::class_<voronoicell_neighbor, boost::noncopyable, py::bases<voronoicell_base>, sptr<voronoicell_neighbor> >(
                //~ "VoronoiCellNeighbors", py::init<>())
        //~ .def("neighbors", vc_neighbors)
    //~ ;
}
