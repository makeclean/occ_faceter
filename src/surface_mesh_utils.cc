#include <iostream>
#include "surface_mesh_utils.hh"

// constructor
surf_utils::surf_utils(std::shared_ptr<moab::Core> MBI, bool with_tagging) {
  moab = MBI;
  tag_data = with_tagging
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
  walk_mesh_and_make_manifolds(starter);
}

eh_t surf_utils::get_neighbour_triangle_by_nodes(const eh_t node1, const eh_t node2, const eh_t triangle_mask) {
  ehVec_t nodes = {node1,node2};
  ehVec_t adjacent_elements;
  moab::ErrorCode rval = moab->get_adjacencies(&nodes[0], nodes.size(), 2, true, 
		  adjacent_elements, moab::Interface::INTERSECT);

  // we should only have two shared by the triangle
  assert(adjacent_elements.size() == 2);

  // return the one that isnt the triangle mask
  if (adjacent_elements[0] == triangle_mask)
    return adjacent_elements[1];
  else 
    return adjacent_elements[0];
}

// get the neighbour elements to my own as a set
ehSet_t surf_utils::get_element_neighbours(const eh_t element){
  ehVec_t adjacent_elements;
  // adjacancies will need to get the 0 d neighbours, then poll the the 3 pairs in turn
  // to build the adjacent elements
  // first get the nodes of the triangle
  ehVec_t nodes; // nodes of the triangle

  // get the nodes on the triangle
  moab::ErrorCode rval = moab->get_adjacencies(&element, 1, 0, true, 
		  nodes, moab::Interface::UNION);

  // now using the three combinations of sides, get the adjacent element
  eh_t neighbour_triangle = get_neighbour_triangle_by_nodes(nodes[0],nodes[1],element);
  adjacent_elements.push_back(neighbour_triangle);
  neighbour_triangle = get_neighbour_triangle_by_nodes(nodes[0],nodes[2],element);
  adjacent_elements.push_back(neighbour_triangle);
  neighbour_triangle = get_neighbour_triangle_by_nodes(nodes[1],nodes[2],element);
  adjacent_elements.push_back(neighbour_triangle);

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
ehSet_t surf_utils::walk_mesh_and_make_manifold(const eh_t starter) {
  ehSet_t manifold = {starter}; // add the start element to the manifold
  ehSet_t neighbours = get_element_neighbours(starter); // get the neighbours
  manifold.insert(neighbours.begin(), neighbours.end()); // add them to the manifold

  // loop over the neighbours
  while (!neighbours.empty()) {
    neighbours = get_elements_neighbours(neighbours, manifold);
    manifold.insert(neighbours.begin(), neighbours.end());
  }
  //. return the manifold
  return manifold  
}

moab::ErrorCode tag_elements_in_set(ehSet_t element_set, int tag_value) {
  
}

// walk the mesh of elements 
void surf_utils::walk_mesh_and_make_manifolds(const eh_t starter) {
  ehSet_t manifold = walk_mesh_and_make_manifold(starter);
  num_manifolds++;
  if (tag_data) tag_elements_in_set(manifold,num_manifolds);

}