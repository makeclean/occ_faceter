#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "brep_faceter.hh"
#include "read_metadata.hh"
#include "BRep_Builder.hxx"
#include "BRepTools.hxx"
#include "moab/GeomQueryTool.hpp"
#include "moab/GeomTopoTool.hpp"

TEST_CASE("Faceting BREP and writing to MOAB", "[faceter]") {

  FacetingTolerance facet_tol(1.e-3);

  MBTool mbtool;
  mbtool.set_faceting_tol_tag(facet_tol.tolerance);

  const char *input_path = "gluedCompSolid.brep";
  const char *metadata_path = "gluedCompSolid_metadata.json";
  const char *output_path = "test_output.h5m";

  TopoDS_Shape shape;
  BRep_Builder builder;

  BRepTools::Read(shape, input_path, builder);

  REQUIRE(!shape.IsNull());

  MaterialsMap materials_map;
  read_metadata(metadata_path, materials_map);

  sew_and_facet(shape, facet_tol, mbtool, materials_map);

  std::vector<moab::EntityHandle> triangles = mbtool.get_entities_by_dimension(0, 2, true);
  CHECK(triangles.size() == 22);

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
  moab::EntityHandle vol2 = gqt.gttool()->entity_by_id(dim, 2);

  int result = 0;
  double xyz[3] = {1.0, 1.0, 1.0};
  ret = gqt.point_in_volume(vol1, xyz, result);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(result == 1);

  ret = gqt.point_in_volume(vol2, xyz, result);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(result == 0);

  xyz[0] = 11;
  result = 0;
  ret = gqt.point_in_volume(vol2, xyz, result);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(result == 1);

  double start[3] = {1.0, 1.0, 1.0};
  double next_surf_dist;
  moab::EntityHandle next_surf;
  double dir[3] = {-1.0, 0.0, 0.0};
  ret = gqt.ray_fire(vol1, start, dir, next_surf, next_surf_dist);
  REQUIRE(ret == moab::MB_SUCCESS);
  CHECK(next_surf_dist == 1.0);
}
