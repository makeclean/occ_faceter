#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

#include <vector>
#include <array>

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
  moab::ErrorCode make_new_volume(moab::EntityHandle &volume);
  moab::ErrorCode add_surface(moab::EntityHandle volume,
			      facet_data facetData);
  moab::ErrorCode add_surface(moab::EntityHandle volume,
			      facet_data facetData,
			      std::vector<edge_data> edgeData);
  void write_geometry(std::string filename);
  
private:
  moab::ErrorCode check_vertex_exists(std::array<double,3> coord, moab::EntityHandle &tVertex);
  moab::ErrorCode make_new_surface(moab::EntityHandle &surface);
  moab::ErrorCode add_facets_to_surface(moab::EntityHandle,
					facet_data facetData);
  moab::ErrorCode add_facets_and_curves_to_surface(moab::EntityHandle,
						    facet_data facetData,
						    std::vector<edge_data> edge_data);
  private:
  moab::Core *mbi = NULL;
  int volID;
  int surfID;
  moab::EntityHandle rootset;
  moab::Tag geometry_dimension_tag, id_tag;
  moab::Tag faceting_tol_tag, geometry_resabs_tag;
  moab::Tag category_tag;
};
