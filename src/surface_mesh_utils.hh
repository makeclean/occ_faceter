#ifndef SURF_UTILS_HH
#define SURF_UTILS_HH 

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
  // go ahead and make our manifolds
  void walk_mesh_and_make_manifolds(const eh_t starter);
  private:
  std::shared_ptr<moab::Core> moab;
};

#endif // SURF_UTILS_HH
