#include "MBTool.hpp"
#include <iostream>
#include <cstring>

#include "moab/GeomTopoTool.hpp"

const char geom_categories[][CATEGORY_TAG_SIZE] = {"Vertex\0",
						   "Curve\0",
						   "Surface\0",
						   "Volume\0",
						   "Group\0"};



static inline void check_moab_rval(moab::ErrorCode ret, const char * what) {
  if (ret == moab::MB_SUCCESS) {
    return;
  }
  throw mberror(ret, what);
}

// default constructor
MBTool::MBTool() {
    mbi = new moab::Core();
    geom_tool = new moab::GeomTopoTool(mbi);

    // new vertex inserter
    vi = new VertexInserter::VertexInserter(mbi,1.e-6); // should pass the
                                        // tolernace by arg  
    for (int i = 0; i < 5; i++) {
      entity_id[i] = 0;
    }
    degenerate_triangle_count = 0;

    // make a new meshset to put stuff in
    moab::ErrorCode rval = mbi->create_meshset(moab::MESHSET_SET, rootset);

    // create tags
    rval = mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                               geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");

    rval = mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                               id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");

    rval = mbi->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                               name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating name_tag");

    rval = mbi->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
                               category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating category_tag");

    rval = mbi->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
                               moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating faceting_tol_tag");

    rval = mbi->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
                               geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating geometry_resabs_tag");

    rval = mbi->tag_get_handle("VOL_ID", 1, moab::MB_TYPE_INTEGER, vol_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get vol_id tag");

    rval = mbi->tag_get_handle("SURF_ID", 1, moab::MB_TYPE_INTEGER, surf_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get surf_id tag");

    rval = mbi->tag_get_handle("MatID", 1, moab::MB_TYPE_INTEGER,
                               mat_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating mat_id_tag");
}

// destructore
MBTool::~MBTool() {
    delete mbi;
    delete geom_tool;
    delete vi;
}

void MBTool::set_faceting_tol_tag(double faceting_tol) {
  moab::EntityHandle set = rootset; // ? *rootset : 0;
  check_moab_rval(mbi->tag_set_data(faceting_tol_tag, &set, 1, &faceting_tol), "set faceting tolerance");
}

void MBTool::set_geometry_tol_tag(double geom_tol) {
  moab::EntityHandle set = rootset; // ? *rootset : 0;
  check_moab_rval(mbi->tag_set_data(geometry_resabs_tag, &set, 1, &geom_tol), "set geometry tol");
}

moab::EntityHandle MBTool::create_entity_set(int dim) {
  // TODO: Check meshset options - Cubit-plugin has different behaviour
  unsigned int options = (dim == 2 || dim == 3) ? moab::MESHSET_ORDERED : moab::MESHSET_SET;

  moab::EntityHandle result;
  check_moab_rval(mbi->create_meshset(options, result), "create entity meshset");

  entity_id[dim]++;
  check_moab_rval(mbi->tag_set_data(id_tag,&result,1,&entity_id[dim]), "create entity set id tag");

  if (dim <= 3) {
    check_moab_rval(mbi->tag_set_data(geometry_dimension_tag,&result,1,&dim), "create entity set dimension tag");
  }

  check_moab_rval(mbi->tag_set_data(category_tag,&result,1,&geom_categories[dim]), "create entity set category tag");
  return result;
}

moab::EntityHandle MBTool::make_new_volume() {
  return create_entity_set(3);
}

moab::EntityHandle MBTool::make_new_surface() {
  return create_entity_set(2);
}

moab::EntityHandle MBTool::make_new_curve() {
  return create_entity_set(1);
}

moab::EntityHandle MBTool::make_new_node() {
  return create_entity_set(0);
}

// add a new group (for materials)
void MBTool::add_group(
  const std::string &name,
  const std::vector<moab::EntityHandle> &entities)
{
  moab::EntityHandle group = create_entity_set(4);

  char namebuf[NAME_TAG_SIZE];
  memset(namebuf, '\0', NAME_TAG_SIZE);
  strncpy(namebuf, name.c_str(), NAME_TAG_SIZE - 1);
  if (name.length() >= (unsigned)NAME_TAG_SIZE) {
    std::cout << "WARNING: group name '" << name.c_str()
              << "' truncated to '" << namebuf << "'" << std::endl;
  }
  check_moab_rval(mbi->tag_set_data(name_tag, &group, 1, namebuf), "set name of added group");
  check_moab_rval(mbi->add_entities(group, entities.data(), entities.size()), "add entities to group");
}

// set a MatID on every triangle (in prep for converting output to .vtk, and viewing in Paraview)
void MBTool::add_mat_ids() {
  // get material groups
  moab::Range material_groups;
  const void *tag_data = &geom_categories[4];
  check_moab_rval(mbi->get_entities_by_type_and_tag(
    0, moab::MBENTITYSET, &category_tag, &tag_data, 1, material_groups),
    "get entities by type and tag");

  for (auto &&group : material_groups) {
    // get groupID
    int groupID;
    check_moab_rval(mbi->tag_get_data(id_tag, &group, 1, &groupID), "get group id");

    // get volumes
    moab::Range volumes;
    check_moab_rval(mbi->get_entities_by_handle(group, volumes), "get group volumes");

    // set MatID on triangles
    for (auto &&vol : volumes) {
      moab::Range surfaces;
      check_moab_rval(mbi->get_child_meshsets(vol, surfaces), "get volume surfaces");

      for (auto &&surface : surfaces) {
        moab::Range tris;
        check_moab_rval(mbi->get_entities_by_type(surface, moab::MBTRI, tris, true), "get surface triangles");

        std::vector<int> matIDs(tris.size(), groupID);
        check_moab_rval(mbi->tag_set_data(mat_id_tag, tris, matIDs.data()), "set triangle material");
      }
    }
  }
}

// add surfaces to volumes, or curves to surfaces
void MBTool::add_child_to_parent(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense)
{
  check_moab_rval(mbi->add_parent_child(volume, surface), "add parent child");
  check_moab_rval(geom_tool->set_sense(surface, volume, sense), "set sense");
}

void MBTool::generate_facet_vertex_map(facet_vertex_map& vertex_map,
                                       const facet_coords& coords)
{
  vertex_map.clear();

  // for each coordinate in the surface make the moab vertex
  int idx = 1; // index start at 1!!
  for ( std::array<double,3> coord : coords ) {
    moab::EntityHandle vert;
    // check for the existence of a vertex, and create a new one if necessary
    moab::ErrorCode rval = vi->insert_vertex(coord,vert);

    moab::EntityHandle set;
    if (rval == moab::MB_ENTITY_NOT_FOUND) {
      // create a meshset for the new vertex
      set = make_new_node();

      moab::Range vertices;
      vertices.insert(vert);
      check_moab_rval(mbi->add_entities(set, vertices), "add vertex");

      vertex_to_set_map[vert] = set;
    } else {
      // find the corresponding meshset
      ent_ent_map::iterator it = vertex_to_set_map.find(vert);
      if (vertex_to_set_map.end() == it) {
        throw std::runtime_error("internal error, vertex not found in map");
      }
      set = it->second;
    }

    facet_vertex_pair v_pair;
    v_pair.vertex = vert;
    v_pair.set = set;
    vertex_map[idx] = v_pair;
    idx += 1;
  }
}

// add facets to surface
void MBTool::add_facets_to_surface(moab::EntityHandle surface,
  const facet_connectivity& connectivity_list, const facet_vertex_map& vertex_map)
{
  moab::Range vertices;
  for (auto const & pair : vertex_map) {
    vertices.insert(pair.second.vertex);
  }
  check_moab_rval(mbi->add_entities(surface, vertices), "add verts to surface");

  moab::Range triangles;
  for ( std::array<int,3> connectivity : connectivity_list) {
    moab::EntityHandle tri;
    moab::EntityHandle connections[3];
    connections[0] = vertex_map.at(connectivity[0]).vertex;
    connections[1] = vertex_map.at(connectivity[1]).vertex;
    connections[2] = vertex_map.at(connectivity[2]).vertex;

    if (connections[2] == connections[1] ||
        connections[1] == connections[0] ||
        connections[2] == connections[0] ) {
      degenerate_triangle_count++;
    } else {
      check_moab_rval(mbi->create_element(moab::MBTRI,connections,3,tri), "create triangle");
      triangles.insert(tri);
    }
  }

  check_moab_rval(mbi->add_entities(surface,triangles), "add triangles to surface");

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
}

// add curves to surface
void MBTool::build_curve(moab::EntityHandle curve,
  edge_data edge_collection_i, const facet_vertex_map& vertex_map)
{
  if (edge_collection_i.connectivity.empty()) {
    std::cout << "Warning: Attempting to build empty curve." << std::endl;
    // throw std::runtime_error("attempting to build empty curve");
  }

  moab::Range edges;
  moab::Range vertices;

  int end_point = 0;
  (edge_collection_i.connectivity.size() - 2  == 0 ) ?
  (end_point = edge_collection_i.connectivity.size() - 1) :
  (end_point = edge_collection_i.connectivity.size() - 2); 

  facet_vertex_pair v_pair = vertex_map.at(edge_collection_i.connectivity.at(0));
  moab::EntityHandle connection[2];
  connection[1] = v_pair.vertex;
  vertices.insert(v_pair.vertex);
  check_moab_rval(mbi->add_parent_child(curve,v_pair.set), "add vertex[0] to curve");

  for ( int j = 0 ; j < end_point ; j++ ) {
    connection[0] = connection[1];
    v_pair = vertex_map.at(edge_collection_i.connectivity.at(j+1));
    connection[1] = v_pair.vertex;
    vertices.insert(v_pair.vertex);

    check_moab_rval(mbi->add_parent_child(curve,v_pair.set), "add vertex[j] to curve");

    moab::EntityHandle edge;
    check_moab_rval(mbi->create_element(moab::MBEDGE, connection, 2, edge), "create edge");
    edges.insert(edge);
  }

  // if curve is closed, remove duplicate vertex
  if (vertices.front() == vertices.back())
    vertices.pop_back();

  // add vertices and edges to curve, and curve to surface
  check_moab_rval(mbi->add_entities(curve,vertices), "add verticies to curve");
  check_moab_rval(mbi->add_entities(curve,edges), "add edges to curve");
}

// write the geometry
void MBTool::write_geometry(const std::string &filename) {
  check_moab_rval(mbi->write_file(filename.c_str()), "write file");

  if (degenerate_triangle_count > 0) {
    std::cout << "Warning: " << degenerate_triangle_count
      << " degenerate triangles have been ignored." << std::endl;
  }
}

std::vector<moab::EntityHandle> MBTool::get_entities_by_dimension(
  const moab::EntityHandle meshset, const int dimension,
  const bool recursive) const
{
  std::vector<moab::EntityHandle> entities;
  check_moab_rval(mbi->get_entities_by_dimension(meshset, dimension, entities, recursive), "get entities by dimension");
  return entities;
}

// add all entities to rootset
void MBTool::gather_ents()
{
  std::cout << "ent counts: " << entity_id[0] << " " << entity_id[1]
  << " " << entity_id[2] << " " << entity_id[3] << " " << entity_id[4] << std::endl;


  moab::Range ents;
  check_moab_rval(mbi->get_entities_by_handle(0, ents), "get entities by handle");
  check_moab_rval(mbi->add_entities(rootset,ents), "add root entities");
}
