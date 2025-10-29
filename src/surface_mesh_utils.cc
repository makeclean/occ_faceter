#include <iostream>
#include "surface_mesh_utils.hh"

// constructor
surf_utils::surf_utils(std::shared_ptr<moab::Core> MBI) {
  moab = MBI;
}

// destructor
surf_utils::~surf_utils(){
}

// setup or check the associated state with 
// volume surface sets and proceed to make 
// n-manifolds by walking the connectivity
void surf_utils::make_manifolds() {
  ehVec_t handles;
  moab::ErrorCode rval = moab->get_entities_by_type(0,moab::MBTRI,handles);
  eh_t starter = handles[0];
  std::cout << handles.size() << std::endl;
  walk_mesh_and_make_manifolds(starter);
}

// get the neighbour elements to my own as a set
ehSet_t surf_utils::get_element_neighbours(const eh_t element){
  std::vector<moab::EntityHandle> adjacent_elements;
  // adjacancies will need to get the 0 d neighbours, then poll the the 3 pairs in turn
  // to build the adjacent elements
  moab::ErrorCode rval = moab->get_adjacencies(&element, 1, 0, true, 
		  adjacent_elements, moab::Interface::UNION);
  std::cout << element << " " << adjacent_elements[0] << " " << adjacent_elements.size() << std::endl;
  std::set<moab::EntityHandle> adjacent_elements_set(adjacent_elements.begin(),
						       adjacent_elements.end());
  return adjacent_elements_set;
}

// get the neighbour elements of vector of elements
ehSet_t surf_utils::get_elements_neighbours(const ehVec_t elements, const ehSet_t exclusions) {
  ehSet_t adjacent_elements;
  for ( moab::EntityHandle element : elements ) {
    ehSet_t neighbours = get_element_neighbours(element);
    adjacent_elements.insert(neighbours.begin(), neighbours.end());
  }
  // eventually we should have visited all elements and at some point we will end up
  // returning an empty set
  if(!exclusions.empty()) {
    std::set<moab::EntityHandle> unique_adjacent_elements;
    std::set_difference(adjacent_elements.begin(), adjacent_elements.end(),
                        exclusions.begin(), exclusions.end(),
                        std::inserter(unique_adjacent_elements, 
                          unique_adjacent_elements.begin()));
    // return the list of elements that are not included in 
    // the exclusions                     
    return unique_adjacent_elements;
  }
  
  return adjacent_elements;
}

// get the neighbour elements of vector of elements
ehSet_t surf_utils::get_elements_neighbours(const ehSet_t elements, const ehSet_t exclusions) {
  ehVec_t elements_vec;
  std::copy(elements.begin(), elements.end(), std::back_inserter(elements_vec));
  return get_elements_neighbours(elements_vec, exclusions);
}

// walk the mesh of elements 
void surf_utils::walk_mesh_and_make_manifolds(const eh_t starter) {
            
  ehSet_t manifold = {starter};
  int delta = 1;
  ehSet_t neighbours = get_element_neighbours(starter);
  manifold.insert(neighbours.begin(), neighbours.end());
  std::cout << "size: " << manifold.size() << std::endl;
  for ( eh_t element : manifold ) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
  // loop over the neighbours
  while (!neighbours.empty()) {
    neighbours = get_elements_neighbours(neighbours, manifold);
    manifold.insert(neighbours.begin(), neighbours.end());
  }

  //
  std::cout << neighbours.size() << " " << manifold.size() << std::endl;
}
