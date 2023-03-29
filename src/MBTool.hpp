#ifndef MBTOOL_HPP
#define MBTOOL_HPP

#include <array>
#include <map>
#include <vector>

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "xyz_to_entity_map.hh"

namespace moab {
class GeomTopoTool;
}

class mberror : public std::runtime_error {
public:
  mberror(moab::ErrorCode error_code, const std::string & what)
    : error_code(error_code), std::runtime_error(what) {}

  moab::ErrorCode error_code;
};

typedef std::vector<moab::EntityHandle> entity_vector;

typedef std::vector<std::array<double, 3>> facet_coords;
typedef std::vector<std::array<int, 3>> facet_connectivity;
typedef std::map<moab::EntityHandle, moab::EntityHandle> ent_ent_map;

struct edge_data {
  std::vector<int> connectivity;
};

class MBTool {
public:
  MBTool();
  ~MBTool();

  void set_scale_factor(double scale_factor);
  void set_faceting_tol_tag(double faceting_tol);
  void set_geometry_tol_tag(double geom_tol);
  moab::EntityHandle make_new_volume();
  moab::EntityHandle make_new_surface();
  moab::EntityHandle make_new_curve();
  moab::EntityHandle make_new_vertex();

  void write_geometry(const std::string &filename);

  moab::EntityHandle find_or_create_vertex(std::array<double,3> point);
  moab::EntityHandle create_triangle(std::array<moab::EntityHandle,3> verticies);
  moab::EntityHandle create_edge(std::array<moab::EntityHandle,2> verticies);
  void add_entities(moab::EntityHandle parent, const std::vector<moab::EntityHandle> &children);

  void note_degenerate_triangle() {
    degenerate_triangle_count += 1;
  }

  void add_node_to_meshset(moab::EntityHandle meshset,
                           std::array<double,3> coord);
  void add_child_to_parent(moab::EntityHandle child,
                           moab::EntityHandle parent, int sense);
  void add_child_to_parent(moab::EntityHandle child,
                           moab::EntityHandle parent);
  void add_group(const std::string &name,
                 const std::vector<moab::EntityHandle> &entities);
  void add_mat_ids();

  std::vector<moab::EntityHandle> get_entities_by_dimension(
      const moab::EntityHandle meshset, const int dimension,
      const bool recursive) const;
  size_t get_number_of_meshsets();
  void gather_ents();

private:
  moab::EntityHandle create_entity_set(int dim);

  moab::Core *mbi;
  coordinates_to_entity_map verticies;

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
  double scale_factor;
};
#endif // MBTOOL_HPP
