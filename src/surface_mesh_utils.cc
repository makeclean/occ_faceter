#include <iostream>
#include <set>
#include <algorithm>
#include "moab/Core.hpp"

// get the neighbour elements to my own as a set
std::set<moab::EntityHandle> get_element_neighbours(const moab::EntityHandle element){
  std::vector<moab::EntityHandle> adjacent_elements;
  moab::ErrorCode rval = moab->get_adjacancies(element, 2, true, &adjacent_elements);
  std::set<moab::EntityHandle> adjacent_elements_set(adjacent_elements.begin(),
						       adjacent_elements.end());
  return adjacent_elements_set;
}

// get the neighbour elements of vector of elements
std::set<moab::EntityHandle> get_elements_neighbours(const std::vector<moab::EntityHandle> elements, const std::set<moab::EntityHandle> exclusions) {
  std::set<moab::EntityHandle> adjacent_elements;
  for ( moab::EntityHandle element : elements ) {
    std::vector<moab::EntityHandle> neighbours = get_element_neighbours(element);
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
std::set<moab::EntityHandle> get_elements_neighbours(const std::set<moab::EntityHandle> elements, const std::set<moab::EntityHandle> exclusions) {
  std::vector<moab::EntityHandle> elements_vec;
  std::copy(elements.begin(), elements.end(), std::back_inserter(elements_vec));
  return get_elements_neighbours(elements_vec, exclusions);
}

// walk the mesh of elements 
void walk_mesh_and_make_manifolds(moab::Core *moab,
				  const moab::EntityHandle starter) {
            
  std::set<moab::EntityHandle> manifold = {starter};
  int delta = 1;
  std::set<moab::EntityHandle> neighbours = get_element_neighbours(starter);
  manifold.insert(neighbours.begin(), neighbours.end());

  // loop over the neighbours
  while (!neighbours.empty()) {
    neighbours = get_elements_neighbours(neighbours, manifold);
    manifold.insert(neighbours.begin(), neighbours.end());
  }

  //
}
