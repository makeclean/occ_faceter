#include "dagmc_topology.hpp"
#include "MBTagConventions.hpp"
#include <iostream>

// constructor
DAGMCTopology::DAGMCTopology(moab::Core *MBI) {
  if (MBI == NULL) {
    mbi = new moab::Core();
    own_moab = true;
  } else {
    mbi = MBI;
    own_moab = false;
  }
}

// destructor
DAGMCTopology::~DAGMCTopology() {
  if(own_moab) delete mbi;
}

// load the file
moab::ErrorCode DAGMCTopology::load_file(const std::string filename) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = mbi->create_meshset(moab::MESHSET_SET, input_set);
  MB_CHK_SET_ERR(rval,"failed to create meshset");
  rval = mbi->load_file(filename.c_str(),&input_set);
  MB_CHK_SET_ERR(rval,"failed to load the file");
  
  if(rval != moab::MB_SUCCESS) return rval;

  return setup_tags();
}

// setup the tags needed for later query
moab::ErrorCode DAGMCTopology::setup_tags() {
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Couldnt get geom dim tag");
  rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Couldnt get id tag");
  rval = mbi->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
			       moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating faceting_tol_tag");
  rval = mbi->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
			     geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating geometry_resabs_tag");
  /*
  rval = mbi->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
			     category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR(rval, "Error creating category_tag");
  */
  return rval;
}

// identify and fix up coincident surfaces
moab::ErrorCode DAGMCTopology::perform_merge() {
  std::vector<merge_pairs_t> coincident;
  moab::ErrorCode rval = identify_coincident(coincident);
  if ( rval != moab::MB_SUCCESS ) return rval;
  if ( coincident.size() == 0 ) std::cout << "no curve pairs found" << std::endl;
  return moab::MB_SUCCESS;

  return merge_coincident(coincident);
}

// identify coincident surfaces
moab::ErrorCode DAGMCTopology::identify_coincident(std::vector<merge_pairs_t> &coincident) {
  // get the list of curves present and compare them for equivlance
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::Range curve_set;

  // maybe use type and tag
  const int dim = 1;
  const void *tag_value[] = {&dim};
  rval = mbi->get_entities_by_type_and_tag(0,moab::MBENTITYSET,&geometry_dimension_tag,
					   tag_value,1,curve_set);

  MB_CHK_SET_ERR(rval,"couldnt get entities by type and tag");
  // first sort by entity set size - i.e. curves of equal length

  std::map<int,moab::Range> curves_and_counts;
  rval = get_curve_lengths(curve_set, curves_and_counts);
  MB_CHK_SET_ERR(rval,"error from get_curve_lengths");
  
  // then compare by vertex handle - these are already shared

  // list of curves that are equivalent - their parents are surfaces that are the identical
  
  return moab::MB_SUCCESS;
}

// count the number of edges in each curve set and store in the map
moab::ErrorCode DAGMCTopology::get_curve_lengths(moab::Range curves,
						 std::map<int,moab::Range> &curves_and_counts) {
  moab::ErrorCode rval = moab::MB_FAILURE;
  std::cout << curves.size() << " curves" << std::endl;
  for ( moab::EntityHandle curve : curves ) {
    // get the id
    int id[1];
    rval = mbi->tag_get_data(id_tag,&(curve),1,id);

    // curve sets
    moab::Range edges;
    edges.clear();
    rval = mbi->get_entities_by_type(curve,moab::MBEDGE,edges);
    int count = edges.size();
    
    if( count == 0 ) {
      std::cout << "curve " << id[0] << " " << curve << " has no edges." << std::endl;
    } else {
      std::cout << "Curve " << id[0] << " " << curve << " has " << edges.size() << " edges." << std::endl;
    }
    moab::Range existing_range = curves_and_counts[count];
    existing_range.merge(edges);

    // add to the range;
    curves_and_counts[count] = existing_range;
  }
  return moab::MB_SUCCESS;
}

// merge surfaces that have been identified as being identical remove one of them
moab::ErrorCode DAGMCTopology::merge_coincident(const std::vector<merge_pairs_t> coincident, const bool del) {
  return moab::MB_SUCCESS;
}
