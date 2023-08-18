#include "vertex_inserter.hh"
#include <iostream>

namespace VertexInserter {

// global vector of hits
std::vector<int> hits;

bool callback(int id) {  
  hits.push_back(id);
  return true;
}

// constructor
VertexInserter::VertexInserter(moab::Core *moab_ptr, double tolerance) {
  mbi = moab_ptr;
  box_tolerance = tolerance;
  count = 0;
}

// insert the given vertex into the moab database - creates a new
// vertex if it doesnt exist - otherwise returns (by arg) the handle
// to the existing one - rval is ENTITY_NOT_FOUND if its new -
// otherwise MB_SUCCESS if found - MB_FAILURE otherwise
moab::ErrorCode VertexInserter::insert_vertex(std::array<double,3> point,
					      moab::EntityHandle &handle) {
  // find any vertices that are 'close' to the point we
  // are trying to insert
  // make a box and keep it 
  Box search = Box(point,box_tolerance);

  // search the tree
  moab::EntityHandle eh = 0;
  moab::ErrorCode rval = search_tree(search,eh);

  if ( rval == moab::MB_ENTITY_NOT_FOUND ) { // vertex doesnt exist
    // make new vertex
    moab::ErrorCode ec = mbi->create_vertex(point.data(),eh);
    MB_CHK_SET_ERR(ec,"Could not create vertex.");
    // set the eh
    search.setHandle(eh); // the curently existing search box belongs to this handle
    // keep it for later
    boxes.push_back(search);
    rtree.Insert(search.min.data(),search.max.data(),count);
    count++;
    //count++;
  }
  handle = eh; /// will be set to 0 under failure conditions
  return rval;
}

// compare the coords to the vertx pointed to, match = true not = false
bool VertexInserter::compare_vertex(const std::array<double,3> coord,
				    const moab::EntityHandle vert) {
  double xyz[3]; // coords of the vertex
  moab::ErrorCode rval = mbi->get_coords(&vert,1,&xyz[0]);

  if ( coord.data()[0] == xyz[0] && coord.data()[1] == xyz[1] && coord.data()[2] == xyz[2] ) {
    return true;
  } 
  return false;
}

// search the tree for the box
moab::ErrorCode VertexInserter::search_tree(const Box search,
					    moab::EntityHandle &hit) {
  //
  moab::EntityHandle rval = moab::MB_FAILURE;
  int nhits = 0;
  // search the tree
  hits.clear();
  nhits = rtree.Search(search.min.data(),search.max.data(),callback);
  
  if ( !nhits ) {
    return moab::MB_ENTITY_NOT_FOUND;
   } else {
    // check the near hits for exactness
    for ( int  i = 0 ; i < nhits ; i++ ) {
      moab::EntityHandle eh = boxes[hits[i]].getHandle();
      if(compare_vertex(search.centre,eh)) {
	hit = eh;
	return moab::MB_SUCCESS;
	hits.clear();
      }
    }
    hits.clear();
    return moab::MB_ENTITY_NOT_FOUND;
  }
  return moab::MB_FAILURE;
}

}