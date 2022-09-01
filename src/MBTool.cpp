#include "MBTool.hpp"
#include <iostream>
#include <cstring>

#include "moab/GeomTopoTool.hpp"
#include "moab/Types.hpp"

const char geom_categories[][CATEGORY_TAG_SIZE] = {"Vertex\0",
						   "Curve\0",
						   "Surface\0",
						   "Volume\0",
						   "Group\0"};


// throw an mberror with message formatted to include filename, linenumber, and failing expression
static __attribute__((noreturn)) void
raise_moab_error(moab::ErrorCode rval, const char *file, const int line, const char *expr) {
  std::ostringstream err;
  err
    << file << ':' << line
    << " got " << moab::ErrorCodeStr[rval]
    << " from: " << expr;
  throw mberror(rval, err.str());
}

// helper to throw an exception fron a moab::ErrorCode
#define CHECK_MOAB_RVAL(expr) \
  do { moab::ErrorCode rval = (expr); \
    if (rval != moab::MB_SUCCESS) raise_moab_error(rval, __FILE__, __LINE__, #expr); \
  } while (0)


// default constructor
MBTool::MBTool() {
  mbi = nullptr;
  geom_tool = nullptr;
  vi = nullptr;

  for (int i = 0; i < 5; i++) {
    entity_id[i] = 0;
  }
  degenerate_triangle_count = 0;

  try {
    mbi = new moab::Core();
    geom_tool = new moab::GeomTopoTool(mbi);

    // new vertex inserter
    vi = new VertexInserter::VertexInserter(mbi,1.e-6); // should pass the tolernace by arg

    // make a new meshset to put stuff in
    CHECK_MOAB_RVAL(mbi->create_meshset(moab::MESHSET_SET, rootset));

    // create tags
    CHECK_MOAB_RVAL(mbi->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                               geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER,
                               id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                               name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                               category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
                               moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE,
                               geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle("VOL_ID", 1, moab::MB_TYPE_INTEGER, vol_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle("SURF_ID", 1, moab::MB_TYPE_INTEGER, surf_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT));

    CHECK_MOAB_RVAL(mbi->tag_get_handle("MatID", 1, moab::MB_TYPE_INTEGER,
                               mat_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT));
  } catch (...) {
    if (vi) delete vi;
    if (geom_tool) delete geom_tool;
    if (mbi) delete mbi;
    throw;
  }
}

// destructor
MBTool::~MBTool() {
  // destroy in reverse order from creation
  delete vi;
  delete geom_tool;
  delete mbi;
}

void MBTool::set_faceting_tol_tag(double faceting_tol) {
  moab::EntityHandle set = rootset; // ? *rootset : 0;
  CHECK_MOAB_RVAL(mbi->tag_set_data(faceting_tol_tag, &set, 1, &faceting_tol));
}

void MBTool::set_geometry_tol_tag(double geom_tol) {
  moab::EntityHandle set = rootset; // ? *rootset : 0;
  CHECK_MOAB_RVAL(mbi->tag_set_data(geometry_resabs_tag, &set, 1, &geom_tol));
}

moab::EntityHandle MBTool::create_entity_set(int dim) {
  // TODO: Check meshset options - Cubit-plugin has different behaviour
  unsigned int options = (dim == 2 || dim == 3) ? moab::MESHSET_ORDERED : moab::MESHSET_SET;

  moab::EntityHandle result;
  CHECK_MOAB_RVAL(mbi->create_meshset(options, result));

  entity_id[dim]++;
  CHECK_MOAB_RVAL(mbi->tag_set_data(id_tag,&result,1,&entity_id[dim]));

  if (dim <= 3) {
    CHECK_MOAB_RVAL(mbi->tag_set_data(geometry_dimension_tag,&result,1,&dim));
  }

  CHECK_MOAB_RVAL(mbi->tag_set_data(category_tag,&result,1,&geom_categories[dim]));
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
  CHECK_MOAB_RVAL(mbi->tag_set_data(name_tag, &group, 1, namebuf));
  CHECK_MOAB_RVAL(mbi->add_entities(group, entities.data(), entities.size()));
}

// set a MatID on every triangle (in prep for converting output to .vtk, and viewing in Paraview)
void MBTool::add_mat_ids() {
  // get material groups
  moab::Range material_groups;
  const void *tag_data = &geom_categories[4];
  CHECK_MOAB_RVAL(mbi->get_entities_by_type_and_tag(
    0, moab::MBENTITYSET, &category_tag, &tag_data, 1, material_groups));

  for (auto &&group : material_groups) {
    // get groupID
    int groupID;
    CHECK_MOAB_RVAL(mbi->tag_get_data(id_tag, &group, 1, &groupID));

    // get volumes
    moab::Range volumes;
    CHECK_MOAB_RVAL(mbi->get_entities_by_handle(group, volumes));

    // set MatID on triangles
    for (auto &&vol : volumes) {
      moab::Range surfaces;
      CHECK_MOAB_RVAL(mbi->get_child_meshsets(vol, surfaces));

      for (auto &&surface : surfaces) {
        moab::Range tris;
        CHECK_MOAB_RVAL(mbi->get_entities_by_type(surface, moab::MBTRI, tris, true));

        std::vector<int> matIDs(tris.size(), groupID);
        CHECK_MOAB_RVAL(mbi->tag_set_data(mat_id_tag, tris, matIDs.data()));
      }
    }
  }
}

// add surfaces to volumes, or curves to surfaces
void MBTool::add_child_to_parent(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense)
{
  CHECK_MOAB_RVAL(mbi->add_parent_child(volume, surface));
  CHECK_MOAB_RVAL(geom_tool->set_sense(surface, volume, sense));
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
      CHECK_MOAB_RVAL(mbi->add_entities(set, vertices));

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
  CHECK_MOAB_RVAL(mbi->add_entities(surface, vertices));

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
      CHECK_MOAB_RVAL(mbi->create_element(moab::MBTRI,connections,3,tri));
      triangles.insert(tri);
    }
  }

  CHECK_MOAB_RVAL(mbi->add_entities(surface,triangles));

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
  CHECK_MOAB_RVAL(mbi->add_parent_child(curve,v_pair.set));

  for ( int j = 0 ; j < end_point ; j++ ) {
    connection[0] = connection[1];
    v_pair = vertex_map.at(edge_collection_i.connectivity.at(j+1));
    connection[1] = v_pair.vertex;
    vertices.insert(v_pair.vertex);

    CHECK_MOAB_RVAL(mbi->add_parent_child(curve,v_pair.set));

    moab::EntityHandle edge;
    CHECK_MOAB_RVAL(mbi->create_element(moab::MBEDGE, connection, 2, edge));
    edges.insert(edge);
  }

  // if curve is closed, remove duplicate vertex
  if (vertices.front() == vertices.back())
    vertices.pop_back();

  // add vertices and edges to curve, and curve to surface
  CHECK_MOAB_RVAL(mbi->add_entities(curve,vertices));
  CHECK_MOAB_RVAL(mbi->add_entities(curve,edges));
}

// write the geometry
void MBTool::write_geometry(const std::string &filename) {
  CHECK_MOAB_RVAL(mbi->write_file(filename.c_str()));

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
  CHECK_MOAB_RVAL(mbi->get_entities_by_dimension(meshset, dimension, entities, recursive));
  return entities;
}

// add all entities to rootset
void MBTool::gather_ents()
{
  std::cout << "ent counts: " << entity_id[0] << " " << entity_id[1]
  << " " << entity_id[2] << " " << entity_id[3] << " " << entity_id[4] << std::endl;

  moab::Range ents;
  CHECK_MOAB_RVAL(mbi->get_entities_by_handle(0, ents));
  CHECK_MOAB_RVAL(mbi->add_entities(rootset,ents));
}
