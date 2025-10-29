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
}
