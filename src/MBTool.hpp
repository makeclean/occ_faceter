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

class MBTool {
public:
  MBTool();
  ~MBTool();

  moab::ErrorCode set_tags();
  moab::ErrorCode add_surface(moab::EntityHandle volume,
			      facet_data facetData);
  moab::ErrorCode add_surface(moab::EntityHandle volume,
			      facet_data facetData,
			      std::vector<edge_data> edgeData);
  moab::ErrorCode add_vertex(std::array<double,3> coords,
			     moab::EntityHandle &vertex,
			     moab::EntityHandle &vertex_set);
  
  void write_geometry(std::string filename);

  void summarise();

  moab::ErrorCode make_new_volume_tags(moab::EntityHandle &volume);

  moab::ErrorCode make_new_surface_tags(moab::EntityHandle &surface);

  moab::ErrorCode make_new_curve_tags(moab::EntityHandle &curve);

  moab::ErrorCode make_new_vertex_tags(moab::EntityHandle &curve);

  moab::ErrorCode add_facets_and_curves_to_surface(moab::EntityHandle,
						    facet_data facetData,
						    std::vector<edge_data> edge_data);

  moab::ErrorCode add_vertex_to_curve(const moab::EntityHandle vertex,
				      const moab::EntityHandle curve);

  moab::ErrorCode add_curve_to_surface(const moab::EntityHandle curve,
					       const moab::EntityHandle surface,
					       const int sense);
    
  moab::ErrorCode add_surface_to_volume(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense);

  moab::ErrorCode add_group(const std::string &name,
                            const std::vector<moab::EntityHandle> &entities);
  moab::ErrorCode add_mat_ids();

  moab::ErrorCode get_entities_by_dimension(const moab::EntityHandle meshset,
                                            const int dimension,
                                            std::vector<moab::EntityHandle> &entities,
                                            const bool recursive) const;
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
  int vertexID;
  moab::EntityHandle rootset;
  moab::Tag geometry_dimension_tag, id_tag;
  moab::Tag faceting_tol_tag, geometry_resabs_tag;
  moab::Tag category_tag;
  moab::Tag vol_id_tag, surf_id_tag; // tags for triangles for plotting
  moab::Tag name_tag;
  moab::Tag mat_id_tag;
  moab::Range existing_vertices;

  std::map<moab::EntityHandle,moab::EntityHandle> vertex2vertexset;

};
#endif // MBTOOL_HPP
