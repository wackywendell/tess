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

py::list vertex_vector_to_list(const std::vector<double>& v){
    py::list l;
    for(uint i=0; i<v.size(); i += 3){
        l.append(py::make_tuple(v[i], v[i+1], v[i+2]));
    }
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
        inline double number_of_edges(){return cell.number_of_edges();}
        
        py::list neighbors(){
            std::vector<int> ns;
            cell.neighbors(ns);
            return std_vector_to_py_list(ns);
        }
        
        py::list vertices(){
            std::vector<double> verts;
            cell.vertices(x,y,z,verts);
            return vertex_vector_to_list(verts);
        }
        
        py::list normals(){
            std::vector<double> verts;
            cell.normals(verts);
            return vertex_vector_to_list(verts);
        }
        
        py::list face_vertices(){
            std::vector<int> f_verts;
            py::list l;
            for(uint i=0; i<f_verts.size(); i += f_verts[i] + 1){
                py::list fs;
                for(uint j=0; j<i; ++j){
                    fs.append(f_verts[i+j+1]);
                }
                l.append(fs);
            }
            return l;
        }
        
        py::list face_areas(){
            std::vector<double> ns;
            cell.face_areas(ns);
            return std_vector_to_py_list(ns);
        }
        
        py::list face_perimeters(){
            std::vector<double> ns;
            cell.face_perimeters(ns);
            return std_vector_to_py_list(ns);
        }
        
        py::list face_orders(){
            std::vector<int> ns;
            cell.face_orders(ns);
            return std_vector_to_py_list(ns);
        }
        
        py::list vertex_orders(){
            std::vector<int> ns;
            cell.vertex_orders(ns);
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

BOOST_PYTHON_MODULE(pyvoro){
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
            .def("number_of_edges", &NeighborCell::number_of_edges)
            .def("vertices", &NeighborCell::vertices, "Cartesian coordinates of the vertices of this cell")
            .def("face_vertices", &NeighborCell::face_vertices, "\
A list of the index of each vertices for each face.\
Each inner list corresponds to one face, with each index in that list corresponding to one vertex.")
            .def("face_areas", &NeighborCell::face_areas)
            .def("face_perimeters", &NeighborCell::face_perimeters)
            .def("normals", &NeighborCell::normals)
            .def("face_orders", &NeighborCell::face_orders, "a list of the number of sides of each face.")
            .def("vertex_orders", &NeighborCell::vertex_orders)
        ;
}
