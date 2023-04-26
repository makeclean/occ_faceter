#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "brep_faceter.hh"
#include "BRepTools.hxx"
#include "moab/GeomQueryTool.hpp"
#include "moab/GeomTopoTool.hpp"
#include "BRepPrimAPI_MakeCylinder.hxx"

TEST_CASE("Faceting cylinder and writing to MOAB", "[cylinder_faceting]") {

  FacetingTolerance facet_tol(1.e-3);

  MBTool mbtool;
  mbtool.set_faceting_tol_tag(facet_tol.tolerance);

  const char *brep_path = "cylinder.brep";
  const char *output_path = "cylinder_test_output.h5m";
  
  float radius = 2;
  float height = 10;
  gp_Ax2 anAxis;
  anAxis.SetLocation(gp_Pnt(0.0, 0, 0.0));
  TopoDS_Shape shape = BRepPrimAPI_MakeCylinder(anAxis, radius, height).Shape();

  BRepTools::Write(shape, brep_path);

  REQUIRE(!shape.IsNull());

  std::vector<std::string> materials_list;
  materials_list.push_back("testmat");

  entity_vector volumes = sew_and_facet2(shape, facet_tol, mbtool);
  add_materials(mbtool, volumes, materials_list);

  std::vector<moab::EntityHandle> triangles = mbtool.get_entities_by_dimension(0, 2, true);

  // This is a check if something has changed - I haven't independantly checked this value:
  CHECK(triangles.size() == 322);

  // 2 vertices + 3 edges + 3 faces + 1 volume + 1 material group + 1 gather set = 11
  // TODO: Add more checks that the meshsets are actually as we expect.
  CHECK(mbtool.get_number_of_meshsets() == 11);

  mbtool.gather_ents();
  mbtool.write_geometry(output_path);

  // pt_in_vol and ray_file tests

  moab::Core mbi;
  moab::ErrorCode ret = mbi.load_file(output_path);
  REQUIRE(ret == moab::MB_SUCCESS);

  moab::GeomTopoTool gtt(&mbi);
  moab::GeomQueryTool gqt(&gtt);
  ret = gqt.initialize();
  REQUIRE(ret == moab::MB_SUCCESS);

  int dim = 3;
  moab::EntityHandle vol1 = gqt.gttool()->entity_by_id(dim, 1);

  int result = 0;
  double xyz[3] = {1.0, 1.0, 1.0};
  ret = gqt.point_in_volume(vol1, xyz, result);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(result == 1);

  double xyz100[3] = {1.0, 1.0, 100.0};
  ret = gqt.point_in_volume(vol1, xyz100, result);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(result == 0);
}
