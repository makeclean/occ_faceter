#ifndef VERTEX_INSERTER_HH
#define VERTEX_INSERTER_HH 1

#include "rtree/RTree.h"
#include "moab/Types.hpp"
#include "moab/Core.hpp"
#include <array>
#include <iostream>

// class to handle the insertion of vertices into moab such
// that they are garenteed to be unique. Each point to be added
// is searched in the tree within the specified tolerance and
// any valid hits are are stored. Each hit is compared for exactness
// if any hit returned is valid it is returned, otherwise the
// vertex is added to the database

// Box struct to assist in the creation of the rtree. Boxes
// store which moab entity handle they surround 
//

namespace VertexInserter {

struct Box {
  // empty constuctor
  Box();
  Box(std::array<double,3> point,
      double tolerance,
      moab::EntityHandle handle = 0) {
    
    centre = point;
    min = { point.data()[0] - tolerance,
	    point.data()[1] - tolerance,
	    point.data()[2] - tolerance};
    max = { point.data()[0] + tolerance,
	    point.data()[1] + tolerance,
	    point.data()[2] + tolerance};
    this->handle = handle;
   }
  
  // array based constructor
  Box(std::array<double,3> centre,
      std::array<double,3> min,
      std::array<double,3> max,
      moab::EntityHandle handle = 0) {
    this->centre = centre;
    this->min = min;
    this->max = max;
    this->handle = handle;
  }
  // simple constructor - calls array version
  Box(double x_0, double y_0, double z_0,
      double x_min, double y_min, double z_min,
      double x_max, double y_max, double z_max,
      moab::EntityHandle handle = 0) {
    std::array<double,3> centre = {x_0,y_0,z_0};
    std::array<double,3> min = {x_min,y_min,z_min};
    std::array<double,3> max = {x_max,y_max,z_max};
    Box(centre,min,max,handle);
  }
  
  // set the handle that this box points to
  void setHandle(moab::EntityHandle handle) {
    this->handle = handle;
  }

  // get the handle this box points to
  moab::EntityHandle getHandle() {
    return handle;
  }

  // print the contents of the box
  void print() {
    std::cout << std::scientific << std::endl;
    std::cout << "boxboxboxboxbox" << std::endl;
    std::cout << "center ";
    std::cout << centre.data()[0] << " ";
    std::cout << centre.data()[1] << " ";
    std::cout << centre.data()[2] << std::endl;
    std::cout << "min " << min.data()[0] << " ";
    std::cout << min.data()[1] << " ";
    std::cout << min.data()[2] << std::endl;
    std::cout << "max " << max.data()[0] << " ";
    std::cout << max.data()[1] << " ";
    std::cout << max.data()[2] << std::endl;
    std::cout << "boxboxboxboxbox" << std::endl;
    return;
  }
  
  public:
    std::array<double,3> centre; /// the midpoint of the box
    std::array<double,3> min;    /// lowest coordinates of the box
    std::array<double,3> max;    /// largest coorinates of the box
  private:
    moab::EntityHandle handle;   /// entity handle of the vertex  
};

class VertexInserter {
  public:
    VertexInserter(moab::Core *mbi, double tolerance = 1.0e-6);

    // insert a vertex into the tree, if rval moab::entity_not_found
    // its a new unique vertex - else moab::entity_found handle points
    // to the existing vertex
    moab::ErrorCode insert_vertex(std::array<double,3> point,
				  moab::EntityHandle &handle);
  
  private:
    // search the rtree for any hits with the specified box
    // returns moab::MB_SUCCESS for found an exisiting entry
    // returns moab::MB_ENTITY_NOT_FOUND for a new entry
    // returns moab::MB_FAILURE for everthing else
    moab::ErrorCode search_tree(const Box search_box,
				moab::EntityHandle &hit);

    // compare the coords stored in coord, with the true
    // coordinates found in vert 
    bool compare_vertex(const std::array<double,3> coord,
			const moab::EntityHandle vert);
  
  private:
    moab::Core *mbi; // class does not own it should not delete
    RTree::RTree<int,double,3,float> rtree; /// the rtree to check against 
    double box_tolerance; /// the insertion search tolerance
    int count; /// the number of boxes in the rtree
    std::vector<Box> boxes; // vector of boxes added
};
}
#endif // VERTEX_INSERTER_HH
