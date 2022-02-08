#ifndef MBTOOL_HPP
#define MBTOOL_HPP

#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <array>

#include "rtree/RTree.h"
#include "vertex_inserter.hh"

namespace moab {
class GeomTopoTool;
}

// convenient return for facets
struct facet_data {
  std::vector<std::array<double,3> > coords;
  std::vector<std::array<int, 3> > connectivity;
};

struct edge_data {
  std::vector<int> connectivity;
};

typedef std::map<int,moab::EntityHandle> facet_vertex_map;

class MBTool {
public:
  MBTool();
  ~MBTool();

  moab::ErrorCode set_tags();
  moab::ErrorCode make_new_volume(moab::EntityHandle &volume);
  void write_geometry(std::string filename);

  void summarise();
  
  moab::ErrorCode make_new_surface(moab::EntityHandle &surface);
  void generate_facet_vertex_map(facet_vertex_map& vertex_map, facet_data facetData);
  moab::ErrorCode add_facets_to_surface(moab::EntityHandle surface,
						    facet_data facetData, const facet_vertex_map& vertex_map);
  moab::ErrorCode add_curves_to_surface(moab::EntityHandle surface,
						    std::vector<edge_data> edge_data, const facet_vertex_map& vertex_map);
  moab::ErrorCode add_surface_to_volume(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense);
  moab::ErrorCode add_group(const std::string &name,
                            const std::vector<moab::EntityHandle> &entities);
  moab::ErrorCode add_mat_ids();

  moab::ErrorCode get_entities_by_dimension(const moab::EntityHandle meshset,
                                            const int dimension,
                                            std::vector<moab::EntityHandle> &entities,
                                            const bool recursive) const;
  moab::ErrorCode make_new_curve(moab::EntityHandle &curve);
private:
  moab::ErrorCode check_vertex_exists(std::array<double,3> coord, moab::EntityHandle &tVertex);
  moab::ErrorCode add_facets_to_surface(moab::EntityHandle,
					facet_data facetData);
  private:
  moab::Core *mbi;
  VertexInserter::VertexInserter *vi;
  moab::GeomTopoTool *geom_tool;
  int groupID;
  int volID;
  int surfID;
  int curveID;
  int degenerate_triangle_count;
  moab::EntityHandle rootset;
  moab::Tag geometry_dimension_tag, id_tag;
  moab::Tag faceting_tol_tag, geometry_resabs_tag;
  moab::Tag category_tag;
  moab::Tag vol_id_tag, surf_id_tag; // tags for triangles for plotting
  moab::Tag name_tag;
  moab::Tag mat_id_tag;
  moab::Range existing_vertices;
};
#endif // MBTOOL_HPP
