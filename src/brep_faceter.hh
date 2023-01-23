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

void sew_and_facet2(TopoDS_Shape &shape, const FacetingTolerance& facet_tol, MBTool &mbtool,
                    std::vector<std::string> &mat_list, std::string single_material = "",
                    bool special_case = false);
void brep_faceter(std::string brep_file, std::string json_file,
                  const FacetingTolerance &facet_tol, std::string h5m_file, bool add_mat_ids);
void read_materials_list(std::string text_file, std::vector<std::string> &mat_list);

#endif // BREP_FACETER_HH
