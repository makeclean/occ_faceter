#include "MBTool.hpp"
#include <iostream>
#include <cstring>

// default constructor
MBTool::MBTool() {
    if(mbi == NULL)
      mbi = new moab::Core();
    volID = 0;
    surfID = 0;
    // make a new meshset to put stuff in
    moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET, rootset);
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
}

// destructor
MBTool::~MBTool() {
    delete mbi;
}

moab::ErrorCode MBTool::set_tags() {
  moab::EntityHandle set = rootset; // ? *rootset : 0;
  moab::ErrorCode rval;
  double faceting_tol = 1.e-4;
  double geom_tol = 1.e-6;
  rval = mbi->tag_set_data(faceting_tol_tag, &set, 1, &faceting_tol);
  rval = mbi->tag_set_data(geometry_resabs_tag, &set, 1, &geom_tol);
  return moab::MB_SUCCESS;
}

// make a new volume meshset 
moab::ErrorCode MBTool::make_new_volume(moab::EntityHandle &volume) {
  volID++;
  // make a new volume set
  moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET,volume);
  // set the id tag
  rval = mbi->tag_set_data(id_tag,&volume,1,&volID);
  // set the dim tag
  int dim = 3;
  rval = mbi->tag_set_data(geometry_dimension_tag,&volume,1,&dim);
  return rval;
}

// add the facet data to a new surface meshset and
// add the parent child links
moab::ErrorCode MBTool::add_surface(moab::EntityHandle volume, facet_data facetData) {
  moab::EntityHandle surface;
  moab::ErrorCode rval = make_new_surface(surface);
  rval = add_facets_to_surface(surface,facetData);
  rval = mbi->add_parent_child(volume,surface);
  return rval;
}

// add surface with edges too
moab::ErrorCode MBTool::add_surface(moab::EntityHandle volume,
				    facet_data facetData,
				    std::vector<edge_data> edges) {
  moab::ErrorCode rval;
  moab::EntityHandle surface;
  rval = make_new_surface(surface);
  std::cout << edges.size() << std::endl;
  rval = add_facets_and_curves_to_surface(surface,facetData,edges);
  rval = mbi->add_parent_child(volume,surface);
}

//  makes a new surface in moab
moab::ErrorCode MBTool::make_new_surface(moab::EntityHandle &surface) {
  surfID++;
  //  moab::EntityHandle surface;
  moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET, surface);
  // set the id tag
  rval = mbi->tag_set_data(id_tag,&surface,1,&surfID);
  // set the dim tag
  int dim = 2;
  rval = mbi->tag_set_data(geometry_dimension_tag,&surface,1,&dim);
  return moab::MB_SUCCESS;  
}

// check for the existence of a vertex - nessessarily linear - unless
// we make a bounding sphere tree of vertices
moab::ErrorCode MBTool::check_vertex_exists(std::array<double,3> coord, moab::EntityHandle &tVertex) {
  // get the vertices
  tVertex = 0;
  moab::Range vertices;
  moab::ErrorCode rval = mbi->get_entities_by_type(0,moab::MBVERTEX,vertices);
  if ( rval != moab::MB_SUCCESS) {
    return rval;
  }
  // check each vertex in turn
  for ( moab::EntityHandle vert : vertices ) {
    double xyz[3];
    // get the coordinate value
    moab::ErrorCode rval = mbi->get_coords(&vert,1,&xyz[0]);
    if ( rval != moab::MB_SUCCESS) {
      return rval;
    }
    if ( coord[0] == xyz[0] && coord[1] == xyz[1] && coord[2] == xyz[2] ) {
      std::cout << "its a me - mario! " << vert << std::endl;
      std::cout << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
      std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
      tVertex = vert;
      // its a match
      return moab::MB_SUCCESS;
    } 
  }
  std::cout << "no match" << std::endl;
  return moab::MB_ENTITY_NOT_FOUND;
}

// add facets to surface
 moab::ErrorCode MBTool::add_facets_to_surface(moab::EntityHandle surface, facet_data facetData) {
   
   std::map<int,moab::EntityHandle> vertex_map;
   // for each coordinate in the surface make the moab vertex
   int idx = 1; // index start at 1!!
   for ( std::array<double,3> coord : facetData.coords ) {
     moab::EntityHandle vert;
     //double coordinate[3];
     //coordinate[
     int matched_index = 0;
     moab::ErrorCode rval = check_vertex_exists(coord, vert);
     if ( rval == moab::MB_ENTITY_NOT_FOUND ) { 
       moab::ErrorCode rval = mbi->create_vertex(coord.data(), vert);
       if(rval != moab::MB_SUCCESS ) {
	 std::cout << "Could not create vertex" << std::endl;
       }
     } else if ( rval == moab::MB_SUCCESS ) {
       std::cout << "vert " << matched_index << " already exists as " << idx  << std::endl;
     }
     vertex_map[idx] = vert;
     idx++;
   }
   
   //now make the triangles
   moab::Range triangles;
   for ( std::array<int,3> connectivity : facetData.connectivity) {
     moab::EntityHandle tri;
     moab::EntityHandle connections[3];
     std::cout << connectivity[0] << " " << connectivity[1] << " " << connectivity[2] << std::endl;
     connections[0] = vertex_map[connectivity[0]];
     connections[1] = vertex_map[connectivity[1]];
     connections[2] = vertex_map[connectivity[2]];
     moab::ErrorCode rval = mbi->create_element(moab::MBTRI,connections,3,tri);
     triangles.insert(tri);
   }
   moab::ErrorCode rval = mbi->add_entities(surface,triangles);
 }

// add facets to surface
moab::ErrorCode MBTool::add_facets_and_curves_to_surface(moab::EntityHandle surface, facet_data facetData, std::vector<edge_data> edge_collection) {
   
   std::map<int,moab::EntityHandle> vertex_map;
   // for each coordinate in the surface make the moab vertex
   int idx = 1; // index start at 1!!
   for ( std::array<double,3> coord : facetData.coords ) {
     moab::EntityHandle vert;
     //double coordinate[3];
     //coordinate[
     moab::ErrorCode rval = check_vertex_exists(coord, vert);
     if ( rval == moab::MB_ENTITY_NOT_FOUND ) { 
       moab::ErrorCode rval = mbi->create_vertex(coord.data(), vert);
     } 
     vertex_map[idx] = vert;
     idx++;
   }
   //now make the triangles
   moab::Range triangles;
   for ( std::array<int,3> connectivity : facetData.connectivity) {
     moab::EntityHandle tri;
     moab::EntityHandle connections[3];
     connections[0] = vertex_map[connectivity[0]];
     connections[1] = vertex_map[connectivity[1]];
     connections[2] = vertex_map[connectivity[2]];
     moab::ErrorCode rval = mbi->create_element(moab::MBTRI,connections,3,tri);
     triangles.insert(tri);
   }
   // now make the edges
   moab::Range edges;
   moab::ErrorCode rval;
   std::cout << edge_collection.size() << std::endl;
   std::cout << edge_collection[0].connectivity.size() << std::endl;
   for ( int i = 0 ; i < edge_collection.size() ; i++ ) {
     moab::EntityHandle curve;
     // set the dimension
     int dim = 1;
     // update the id
     int id = 1 + i;
     // set the name of the meshet
     const char* name = "CURVE\0";
     rval = mbi->create_meshset(moab::MESHSET_SET, curve);
     rval = mbi->tag_set_data(geometry_dimension_tag,&curve,1,&dim);
     rval = mbi->tag_set_data(category_tag, &curve, 1, &name);
     for ( int j = 0 ; j < edge_collection[i].connectivity.size() - 2 ; j++ ) {
       moab::EntityHandle h;
       moab::EntityHandle connection[2];

       connection[0] = vertex_map[edge_collection[i].connectivity[j]];
       connection[1] = vertex_map[edge_collection[i].connectivity[j+1]];
       rval = mbi->add_parent_child(curve,connection[0]);
       rval = mbi->add_parent_child(curve,connection[1]);
       // create the edge type
       rval = mbi->create_element(moab::MBEDGE, connection, 2, h);
       edges.insert(h);
     }
     // add edges to curve
     rval = mbi->add_entities(curve,edges);
     rval = mbi->add_parent_child(surface,curve);
   }
   // tag the edges appropriately 
   // add the edges to the surface set
   rval = mbi->add_entities(surface,edges); 
   rval = mbi->add_entities(surface,triangles);
 }


  // write the geometry
void MBTool::write_geometry(std::string filename) {
  moab::ErrorCode rval = mbi->write_file(filename.c_str());
  return;
}

