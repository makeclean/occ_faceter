#ifndef BREP_FACETER_HH
#define BREP_FACETER_HH

#include <string>

#ifndef OCE
  #include "TopoDS_Shape.hxx"
#else
  #include "TopoDS_Shape.lxx"
#endif

#include "MBTool.hpp"
#include "read_metadata.hh"

struct FacetingTolerance {
  float tolerance;
  bool is_relative;

  FacetingTolerance(float tol, bool is_absolute = false)
    : tolerance(tol), is_relative(!is_absolute) {}
};

void sew_and_facet(TopoDS_Shape &shape, const FacetingTolerance &facet_tol, MBTool &mbtool,
                   MaterialsMap &mat_map, std::string single_material = "");
void brep_faceter(std::string brep_file, std::string json_file,
                  const FacetingTolerance &facet_tol, std::string h5m_file, bool add_mat_ids);

#endif // BREP_FACETER_HH
