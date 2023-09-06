#ifndef BREP_FACETER_HH
#define BREP_FACETER_HH

#include <string>

#include "TopoDS_Shape.hxx"
#include "MBTool.hpp"

struct FacetingTolerance {
  float tolerance;
  bool is_relative;

  FacetingTolerance(float tol, bool is_absolute = false)
      : tolerance(tol), is_relative(!is_absolute) {}
};

entity_vector sew_and_facet2(TopoDS_Shape &shape, const FacetingTolerance &facet_tol, MBTool &mbtool);

void read_materials_list(std::string text_file, std::vector<std::string> &mat_list);
void add_materials(MBTool &mbtool, const entity_vector &volumes, const std::vector<std::string> &mat_list);
void add_single_material(MBTool &mbtool, const entity_vector &volumes, const std::string &single_material);

void brep_faceter(std::string brep_file, std::string json_file,
                  const FacetingTolerance &facet_tol, std::string h5m_file,
                  bool add_mat_ids, double scale_factor);

#endif // BREP_FACETER_HH
