#include "dagmc_topology.hpp"

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
  rval = mbi->load_file(filename.c_str(),input_set);
  MB_CHK_SET_ERR(rval,"failed to load the file");
  
  if(rval != moab::MB_SUCCESS) return;

  return setup_tags();
}

// setup the tags needed for later query
moab::ErrorCode DAGMCTopolgy::setup_tags() {
  rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");
  rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");
  rval = mbi->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
			       moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Error creating faceting_tol_tag");
  rval = mbi->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
			     geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Error creating geometry_resabs_tag");
  rval = mbi->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
			     category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Error creating category_tag");
  return rval;
}
