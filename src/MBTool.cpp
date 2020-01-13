#include "MBTool.hpp"
#include <iostream>
#include <cstring>

#include "moab/GeomTopoTool.hpp"

const char geom_categories[][CATEGORY_TAG_SIZE] = {"Vertex\0",
						   "Curve\0",
						   "Surface\0",
						   "Volume\0",
						   "Group\0"};


// default constructor
MBTool::MBTool() {
    mbi = new moab::Core();
    geom_tool = new moab::GeomTopoTool(mbi);

    // new vertex inserter
    vi = new VertexInserter::VertexInserter(mbi,1.e-6); // should pass the
                                        // tolernace by arg  
    
    volID = 0;
    surfID = 0;
    curveID = 0;
    existing_vertices.clear();
    // make a new meshset to put stuff in
    moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET, rootset);
    rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");
    rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");

    rval = mbi->tag_get_handle("VOL_ID", 1, moab::MB_TYPE_INTEGER, vol_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get vol_id tag");

    rval = mbi->tag_get_handle("SURF_ID", 1, moab::MB_TYPE_INTEGER, surf_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get surf_id tag");
    
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
    delete geom_tool;
    delete vi;
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

  // std::cout << "Created new volume " << volID << std::endl;
  
  // make a new volume set
  moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_ORDERED,volume);
  // set the id tag
  rval = mbi->tag_set_data(id_tag,&volume,1,&volID);
  // set the dim tag
  int dim = 3;
  rval = mbi->tag_set_data(geometry_dimension_tag,&volume,1,&dim);
  rval = mbi->tag_set_data(category_tag,&volume,1,&geom_categories[dim]);
  return rval;
}

// add surface with edges too
moab::ErrorCode MBTool::add_surface(moab::EntityHandle volume,
				    facet_data facetData,
				    std::vector<edge_data> edges) {
  moab::ErrorCode rval;
  moab::EntityHandle surface;
  rval = make_new_surface(surface);
  rval = add_facets_and_curves_to_surface(surface,facetData,edges);
  rval = mbi->add_parent_child(volume,surface);
  return rval;
}

moab::ErrorCode MBTool::add_surface_to_volume(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense) {
  moab::ErrorCode rval;
  rval = mbi->add_parent_child(volume, surface);
  MB_CHK_ERR(rval);
  rval = geom_tool->set_sense(surface, volume, sense);
  return rval;
}

//  makes a new surface in moab
moab::ErrorCode MBTool::make_new_surface(moab::EntityHandle &surface) {
  surfID++;
  //  moab::EntityHandle surface;
  // std::cout << "Created new surface " << surfID << std::endl;
  
  moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_ORDERED, surface);
  // set the id tag
  rval = mbi->tag_set_data(id_tag,&surface,1,&surfID);
  // set the dim tag
  int dim = 2;
  rval = mbi->tag_set_data(geometry_dimension_tag,&surface,1,&dim);
  rval = mbi->tag_set_data(category_tag,&surface,1,&geom_categories[dim]);
  return moab::MB_SUCCESS;  
}

//  makes a new surface in moab
moab::ErrorCode MBTool::make_new_curve(moab::EntityHandle &curve) {
  curveID++;
  //  moab::EntityHandle surface;
  moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET, curve);
  // set the id tag
  rval = mbi->tag_set_data(id_tag,&curve,1,&curveID);
  // set the dim tag
  int dim = 1;
  rval = mbi->tag_set_data(geometry_dimension_tag,&curve,1,&dim);
  // set the name of the meshet
  rval = mbi->tag_set_data(category_tag, &curve, 1, &geom_categories[dim]);
  return moab::MB_SUCCESS;  
}

// check for the existence of a vertex, if it exists return the eh otherwise
// create a new one 
moab::ErrorCode MBTool::check_vertex_exists(std::array<double,3> coord,
					    moab::EntityHandle &tVertex) {
  // get the vertices
  tVertex = 0;

  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = vi->insert_vertex(coord,tVertex);
  return rval;
}

// add facets to surface
moab::ErrorCode MBTool::add_facets_and_curves_to_surface(moab::EntityHandle surface, facet_data facetData, std::vector<edge_data> edge_collection) {
   
   std::map<int,moab::EntityHandle> vertex_map;
   // for each coordinate in the surface make the moab vertex
   int idx = 1; // index start at 1!!
   for ( std::array<double,3> coord : facetData.coords ) {
     moab::EntityHandle vert;
     moab::ErrorCode rval = check_vertex_exists(coord, vert);
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

     if ( connections[2] == connections[1] || 
          connections[1] == connections[0] || 
          connections[2] == connections[0] ) {
       std::cout << "degenerate triangle not created" << std::endl;
     } else {
       moab::ErrorCode rval = mbi->create_element(moab::MBTRI,connections,3,tri);
       triangles.insert(tri);
     }
   }
   // now make the edges
   moab::Range edges;
   moab::ErrorCode rval;

   for ( int i = 0 ; i < edge_collection.size() ; i++ ) {
     
     moab::EntityHandle curve;
     // construct an entity set for the curve
     rval = make_new_curve(curve);
     edges.clear();
     
     int end_point = 0;
     (edge_collection[i].connectivity.size() - 2  == 0 ) ?
     (end_point = edge_collection[i].connectivity.size() - 1) :
     (end_point = edge_collection[i].connectivity.size() - 2); 

     for ( int j = 0 ; j < end_point ; j++ ) {
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
   // add the edges to the surface set
   rval = mbi->add_entities(surface,edges); 
   rval = mbi->add_entities(surface,triangles);
   // tag the triangles with

   // tag the triangles with their volume id 
   int vol_value = volID;
   const void *vol_ptr = &vol_value;
   // tag the triangles with their surface id
   int surf_value = surfID;
   const void *surf_ptr = &surf_value;   

   /*
   for ( moab::EntityHandle triangle : triangles ) {
     rval = mbi->tag_set_data(vol_id_tag, &triangle,1, vol_ptr);
     rval = mbi->tag_set_data(surf_id_tag, &triangle,1, surf_ptr);
   }
   */
   return rval;
}


// write the geometry
void MBTool::write_geometry(std::string filename) {
  moab::ErrorCode rval = mbi->write_file(filename.c_str());
  return;
}

void summarise() {
  moab::ErrorCode rval = moab::MB_FAILURE;
  //  rval = mbi->get_entities_by_type_and_tag(0,moab::MBENTITYSET,
}

moab::ErrorCode MBTool::get_entities_by_dimension(const moab::EntityHandle meshset,
                                                  const int dimension,
                                                  std::vector<moab::EntityHandle> &entities,
                                                  const bool recursive) const {
  return mbi->get_entities_by_dimension(meshset, dimension, entities, recursive);
}
