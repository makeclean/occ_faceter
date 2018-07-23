#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "vertex_inserter.hh"
#include "moab/Types.hpp"

TEST_CASE("test insertion", "[insert_vertex]") {
  VertexInserter *vi = new VertexInserter(1.e-6);
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
  REQUIRE(rval == moab::MB_SUCCESS); // entity already exists
  rval = vi->insert_vertex({0,0,1.-1.e-7},eh); // entity 
  REQUIRE(rval == moab::MB_SUCCESS); // entity doesnt exist
  delete vi;
}

