#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "vertex_inserter.hh"
#include "moab/Types.hpp"
#include "moab/Core.hpp"
#include <random>

TEST_CASE("test insertion", "[insert_vertex]") {
  moab::Core *mbi = new moab::Core();
  VertexInserter::VertexInserter *vi = new VertexInserter::VertexInserter(mbi,1.e-6);
  moab::EntityHandle eh = 0;
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = vi->insert_vertex({0,0,0},eh);
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND);
  rval = vi->insert_vertex({0,1,0},eh);
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND);
  rval = vi->insert_vertex({0,0,1},eh);
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND);
  rval = vi->insert_vertex({0,0,0},eh); 
  REQUIRE(rval == moab::MB_SUCCESS); // entity already exists
  rval = vi->insert_vertex({0,1,0},eh);
  REQUIRE(rval == moab::MB_SUCCESS); // entity already exists
  rval = vi->insert_vertex({0,0,1},eh);
  REQUIRE(rval == moab::MB_SUCCESS);
  rval = vi->insert_vertex({0,0,1.-1.e-1},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-2},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-3},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-4},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-5},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-6},eh); // entity 
  REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND); // entity already exists
  delete vi;
  delete mbi;
}

TEST_CASE("random insertion", "[random_insert]") {
  moab::Core *mbi = new moab::Core();
  VertexInserter::VertexInserter *vi = new VertexInserter::VertexInserter(mbi,1.e-6);
  int num = 100000;

  std::random_device rd;
  std::mt19937_64 lcg(rd());
  lcg.seed(1);
  std::uniform_real_distribution<double> uniform_dist(-10.,10.);

  moab::EntityHandle eh;
  moab::ErrorCode rval;
  for ( int i = 0 ; i < num ; i++ ) {
    double x = uniform_dist(lcg);
    double y = uniform_dist(lcg);
    double z = uniform_dist(lcg);
    rval = vi->insert_vertex({x,y,z},eh);
    REQUIRE(rval == moab::MB_ENTITY_NOT_FOUND);
    REQUIRE(rval != moab::MB_FAILURE);
  }
}
