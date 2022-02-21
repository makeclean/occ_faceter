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
    groupID = 0;
    volID = 0;
    surfID = 0;
    curveID = 0;
    degenerate_triangle_count = 0;
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

    rval = mbi->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                               name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating name_tag");

    rval = mbi->tag_get_handle("MatID", 1, moab::MB_TYPE_INTEGER,
                               mat_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating mat_id_tag");
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

// add a new group (for materials)
moab::ErrorCode MBTool::add_group(const std::string &name,
                                  const std::vector<moab::EntityHandle> &entities) {
  groupID++;
  moab::ErrorCode rval;
  moab::EntityHandle group;
  rval = mbi->create_meshset(moab::MESHSET_SET, group);
  if (moab::MB_SUCCESS != rval) return rval;

  rval = mbi->tag_set_data(id_tag, &group, 1, &groupID);
  if (moab::MB_SUCCESS != rval) return rval;

  rval = mbi->tag_set_data(category_tag, &group, 1, &geom_categories[4]);
  if (moab::MB_SUCCESS != rval) return rval;

  char namebuf[NAME_TAG_SIZE];
  memset(namebuf, '\0', NAME_TAG_SIZE);
  strncpy(namebuf, name.c_str(), NAME_TAG_SIZE - 1);
  if (name.length() >= (unsigned)NAME_TAG_SIZE) {
    std::cout << "WARNING: group name '" << name.c_str()
              << "' truncated to '" << namebuf << "'" << std::endl;
  }
  rval = mbi->tag_set_data(name_tag, &group, 1, namebuf);
  if (moab::MB_SUCCESS != rval) return rval;

  rval = mbi->add_entities(group, entities.data(), entities.size());
  return rval;
}

// set a MatID on every triangle (in prep for converting output to .vtk, and viewing in Paraview)
moab::ErrorCode MBTool::add_mat_ids() {
  moab::ErrorCode rval;

  // get material groups
  moab::Range material_groups;
  const void *tag_data = &geom_categories[4];
  rval = mbi->get_entities_by_type_and_tag(0, moab::MBENTITYSET, &category_tag,
                                           &tag_data, 1, material_groups);
  if (moab::MB_SUCCESS != rval) return rval;

  for (auto &&group : material_groups) {
    // get groupID
    int groupID;
    rval = mbi->tag_get_data(id_tag, &group, 1, &groupID);
    if (moab::MB_SUCCESS != rval) return rval;

    // get volumes
    moab::Range volumes;
    rval = mbi->get_entities_by_handle(group, volumes);
    if (moab::MB_SUCCESS != rval) return rval;

    // set MatID on triangles
    for (auto &&vol : volumes) {
      moab::Range surfaces;
      rval = mbi->get_child_meshsets(vol, surfaces);
      if (moab::MB_SUCCESS != rval) return rval;

      for (auto &&surface : surfaces) {
        moab::Range tris;
        rval = mbi->get_entities_by_type(surface, moab::MBTRI, tris, true);
        if (moab::MB_SUCCESS != rval) return rval;

        std::vector<int> matIDs(tris.size(), groupID);
        rval = mbi->tag_set_data(mat_id_tag, tris, matIDs.data());
        if (moab::MB_SUCCESS != rval) return rval;
      }
    }
  }
  return rval;
}

// add surfaces to volumes, or curves to surfaces
moab::ErrorCode MBTool::add_child_to_parent(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense) {
  moab::ErrorCode rval;
  rval = mbi->add_parent_child(volume, surface);
  MB_CHK_ERR(rval);
  return geom_tool->set_sense(surface, volume, sense);
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

//  makes a new curve in moab
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

void MBTool::generate_facet_vertex_map(facet_vertex_map& vertex_map, facet_data facetData) {
  vertex_map.clear();

  // for each coordinate in the surface make the moab vertex
  int idx = 1; // index start at 1!!
  for ( std::array<double,3> coord : facetData.coords ) {
    moab::EntityHandle vert;
    moab::ErrorCode rval = check_vertex_exists(coord, vert);
    vertex_map[idx] = vert;
    idx++;
  }
}

// add facets to surface
moab::ErrorCode MBTool::add_facets_to_surface(moab::EntityHandle surface,
  facet_data facetData, const facet_vertex_map& vertex_map) {

  moab::Range vertices;
  for (auto const & pair : vertex_map) {
    vertices.insert(pair.second);
  }
  moab::ErrorCode rval = mbi->add_entities(surface, vertices);

  moab::Range triangles;
  for ( std::array<int,3> connectivity : facetData.connectivity) {
    moab::EntityHandle tri;
    moab::EntityHandle connections[3];
    connections[0] = vertex_map.at(connectivity[0]);
    connections[1] = vertex_map.at(connectivity[1]);
    connections[2] = vertex_map.at(connectivity[2]);

    if ( connections[2] == connections[1] ||
        connections[1] == connections[0] ||
        connections[2] == connections[0] ) {
      degenerate_triangle_count++;
    } else {
      moab::ErrorCode rval = mbi->create_element(moab::MBTRI,connections,3,tri);
      triangles.insert(tri);
    }
  }

  rval = mbi->add_entities(surface,triangles);

  /*
  // tag the triangles with their volume id
  int vol_value = volID;
  const void *vol_ptr = &vol_value;
  // tag the triangles with their surface id
  int surf_value = surfID;
  const void *surf_ptr = &surf_value;

  for ( moab::EntityHandle triangle : triangles ) {
    rval = mbi->tag_set_data(vol_id_tag, &triangle,1, vol_ptr);
    rval = mbi->tag_set_data(surf_id_tag, &triangle,1, surf_ptr);
  }
  */
  return rval;
}

// add curves to surface
moab::ErrorCode MBTool::build_curve(moab::EntityHandle curve,
  edge_data edge_collection_i, const facet_vertex_map& vertex_map) {
  
  moab::ErrorCode rval;

  moab::Range edges;
  moab::Range vertices;

  int end_point = 0;
  (edge_collection_i.connectivity.size() - 2  == 0 ) ?
  (end_point = edge_collection_i.connectivity.size() - 1) :
  (end_point = edge_collection_i.connectivity.size() - 2); 

  moab::EntityHandle connection[2];
  connection[1] = vertex_map.at(edge_collection_i.connectivity[0]);
  vertices.insert(connection[1]);
  rval = mbi->add_parent_child(curve,connection[1]);

  for ( int j = 0 ; j < end_point ; j++ ) {

    connection[0] = connection[1];
    connection[1] = vertex_map.at(edge_collection_i.connectivity[j+1]);
    vertices.insert(connection[1]);
    rval = mbi->add_parent_child(curve,connection[1]);

    moab::EntityHandle h;
    rval = mbi->create_element(moab::MBEDGE, connection, 2, h);
    edges.insert(h);
  }

  // if curve is closed, remove duplicate vertex
  if (vertices.front() == vertices.back())
    vertices.pop_back();

  // add vertices and edges to curve, and curve to surface
  rval = mbi->add_entities(curve,vertices);
  rval = mbi->add_entities(curve,edges);
  return rval;
}

// write the geometry
void MBTool::write_geometry(std::string filename) {
  moab::ErrorCode rval = mbi->write_file(filename.c_str());

  if (degenerate_triangle_count > 0) {
    std::cout << "Warning: " << degenerate_triangle_count
      << " degenerate triangles have been ignored." << std::endl;
  }
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
