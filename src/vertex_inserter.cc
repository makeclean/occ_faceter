#include "vertex_inserter.hh"
#include <iostream>

// global vector of hits
std::vector<int> hits;

bool callback(int id) {  
  hits.push_back(id);
  return true;
}

// constructor
VertexInserter::VertexInserter(double tolerance) {
  box_tolerance = tolerance;
  count = 0;
}

// desctructor
VertexInserter::~VertexInserter() {
  boxes.clear();
  rtree.RemoveAll();
}

// insert the given vertex into the moab database
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
    rtree.Insert(search.min.data(),search.max.data(),count);
    // set the eh
    eh = count;
    search.setHandle(eh);
    count++;
    // keep it for later
    boxes.push_back(search);
  } else if (rval == moab::MB_SUCCESS) { // vertex already exists
    handle = eh;
  } 
  return rval;
}

// search the tree for the box
moab::ErrorCode VertexInserter::search_tree(const Box search,
					    moab::EntityHandle &hit) {
  // 
  moab::EntityHandle rval = moab::MB_FAILURE;
  int nhits = 0;
  // search the tree
  nhits = rtree.Search(search.min.data(),search.max.data(),callback);
  
  if ( !nhits ) {
    std::cout << "no hit found" << std::endl;
    return moab::MB_ENTITY_NOT_FOUND;
   } else {
    // check the near hits for exactness
    std::cout << nhits << " hits!" << std::endl;
    for ( int  i = 0 ; i < nhits ; i++ ) {
      moab::EntityHandle eh = boxes[hits[i]].getHandle();
      std::cout << "hit box " << hits[i] << " with handle " << eh << std::endl;
    }
    hits.clear();
    return moab::MB_SUCCESS;
  }
  return moab::MB_FAILURE;
}

