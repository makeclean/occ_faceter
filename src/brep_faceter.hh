#ifndef BREP_FACETER_HH
#define BREP_FACETER_HH

#include <string>

#include "TopoDS_Shape.hxx"
#include "MBTool.hpp"
#include "read_metadata.hh"

void sew_and_facet(TopoDS_Shape &shape, float facet_tol, MBTool &mbtool,
                   MaterialsMap &mat_map, std::string single_material = "");
void brep_faceter(std::string brep_file, std::string json_file,
                  float facet_tol, std::string h5m_file);

#endif // BREP_FACETER_HH
