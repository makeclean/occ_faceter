#ifndef SURF_UTILS_HH
#define SURF_UTILS_HH 

#include <cassert>
#include <iostream>
#include <memory>
#include <set>
#include <vector>
#include "moab/Core.hpp"

typedef moab::EntityHandle eh_t;
typedef std::set<moab::EntityHandle> ehSet_t;
typedef std::vector<moab::EntityHandle> ehVec_t;

class surf_utils {
  public:
   // constructor
   surf_utils(std::shared_ptr<moab::Core> MBI);
   // destructor
   ~surf_utils();
   // walk the mesh and make n-manifolds
   void make_manifolds(); 
  private:
  // given two nodes on a triangle, and a mask triangle to ignore, return
  // the other triangle that is shared by nodes 1 and 2
  eh_t get_neighbour_triangle_by_nodes(const eh_t node1, const eh_t node2, 
    const eh_t triangle_mask);

  // given a singular element, find its neighbour elements by adjacancy 
  ehSet_t get_element_neighbours(const eh_t element);

  // given a vector of elements, and a potential exclusion 
  // set (mask) return a set of elements
  ehSet_t get_elements_neighbours(const ehVec_t elements, 
  		                  const ehSet_t exclusions);
  // given a set of elemens, and a potential exclusion
  // set (mask) return a set of elements
  ehSet_t get_elements_neighbours(const ehSet_t elements, 
		  		  const ehSet_t exclusions);

  // assuming a setup of input state that is required
  // go ahead and find all the manifolds
  void walk_mesh_and_make_manifolds(const eh_t starter);

  // assuming a setup of input state that is required
  // go ahead and make a single manifold given a starter
  // element
  void walk_mesh_and_make_manifolds(const eh_t starter);

  private:
  std::shared_ptr<moab::Core> moab; /// the MOAB core instance
  bool tag_data; /// tag the moab set with various useful output data 
  int num_manifolds; /// the number of manifolds discovered
};

#endif // SURF_UTILS_HH
