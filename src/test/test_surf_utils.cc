#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_NO_POSIX_SIGNALS
#include "catch.hpp"
#include "surface_mesh_utils.hh"

TEST_CASE("Test of surface utils for making manifolds", "[surf_utils]") {

  // pt_in_vol and ray_file tests

  std::shared_ptr<moab::Core> mbi = std::make_shared<moab::Core>();
  moab::ErrorCode rval = mbi->load_file("test.h5m");
  REQUIRE(rval == moab::MB_SUCCESS);

  // make a surf_utils instance
  std::shared_ptr<surf_utils> meshUtils = std::make_shared<surf_utils>(mbi,true);
  // 
  meshUtils->make_manifolds();

  // number of elements tagged with manifold_id should equal
  // the number of triangles in the instance
  ehVec_t triangles;
  // get the triangles
  rval = mbi->get_entities_by_type(0,moab::MBTRI,triangles);
  //ehVec_t tagged_triangles;
  moab::Range tagged_triangles;
  moab::Tag tag;
  // get the tag handle
  rval = mbi->tag_get_handle("MANIFOLD_ID", 1, moab::MB_TYPE_INTEGER, tag,
                               moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT); 

  const int tag_value = 1;
  const void *value[1] = {&tag_value};
  // get the triangles that are tagged 
  // with a manfold id of 1
  rval = mbi->get_entities_by_type_and_tag(0,moab::MBTRI,
		  &tag, value, 1, tagged_triangles); 

  REQUIRE(triangles.size() == tagged_triangles.size());

  rval = mbi->write_file("mesh.h5m"); 
}
