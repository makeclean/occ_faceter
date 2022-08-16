#ifndef MBTOOL_HPP
#define MBTOOL_HPP

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <array>
#include <map>

#include "rtree/RTree.h"
#include "vertex_inserter.hh"

namespace moab {
class GeomTopoTool;
}

class mberror : public std::exception {
public:
  mberror(moab::ErrorCode error_code) : err(error_code) {}

  moab::ErrorCode err;
};

typedef std::vector<std::array<double,3>> facet_coords;
typedef std::vector<std::array<int, 3>> facet_connectivity;
typedef std::map<moab::EntityHandle, moab::EntityHandle> ent_ent_map;

// convenient return for facets
struct facet_data {
  facet_coords coords;
  facet_connectivity connectivity;
};

struct edge_data {
  std::vector<int> connectivity;
};

struct facet_vertex_pair {
  moab::EntityHandle vertex;
  moab::EntityHandle set;
};

typedef std::map<int,facet_vertex_pair> facet_vertex_map;

class MBTool {
public:
  MBTool();
  ~MBTool();

  void set_faceting_tol_tag(double faceting_tol);
  void set_geometry_tol_tag(double geom_tol);
  moab::EntityHandle make_new_volume();
  moab::EntityHandle make_new_surface();
  moab::EntityHandle make_new_curve();
  moab::EntityHandle make_new_node();

  void write_geometry(const std::string &filename);

  moab::ErrorCode generate_facet_vertex_map(facet_vertex_map& vertex_map,
                                            const facet_coords& coords);
  moab::ErrorCode add_facets_to_surface(moab::EntityHandle surface,
                                        const facet_connectivity& connectivity_list,
                                        const facet_vertex_map& vertex_map);
  moab::ErrorCode build_curve(moab::EntityHandle curve, edge_data edge,
                              const facet_vertex_map& vertex_map);
  moab::ErrorCode add_child_to_parent(moab::EntityHandle child,
          moab::EntityHandle parent, int sense);
  void add_group(const std::string &name,
                 const std::vector<moab::EntityHandle> &entities);
  moab::ErrorCode add_mat_ids();

  moab::ErrorCode get_entities_by_dimension(const moab::EntityHandle meshset,
                                            const int dimension,
                                            std::vector<moab::EntityHandle> &entities,
                                            const bool recursive) const;
  moab::ErrorCode gather_ents();

private:
  moab::EntityHandle create_entity_set(int dim);

  moab::Core *mbi;
  VertexInserter::VertexInserter *vi;
  ent_ent_map vertex_to_set_map;

  moab::GeomTopoTool *geom_tool;
  int entity_id[5]; // group, volume, surface, curve IDs (indexed by dim)
  int degenerate_triangle_count;
  moab::EntityHandle rootset;
  moab::Tag geometry_dimension_tag, id_tag;
  moab::Tag faceting_tol_tag, geometry_resabs_tag;
  moab::Tag category_tag;
  moab::Tag vol_id_tag, surf_id_tag; // tags for triangles for plotting
  moab::Tag name_tag;
  moab::Tag mat_id_tag;
};
#endif // MBTOOL_HPP
