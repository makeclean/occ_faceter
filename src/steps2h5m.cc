#include "steps2h5m.hh"

#include <iostream>
#include <fstream>

#include "STEPControl_Reader.hxx"
#include "brep_faceter.hh"
#include "MBTool.hpp"
#include "step2breps.hh"
#include "UniqueId/json.hpp"

using json = nlohmann::json;

void steps2h5m(std::string input_file, const FacetingTolerance &facet_tol, std::string h5m_file) {

  json j = json::parse(std::ifstream(input_file));

  MBTool mbtool;
  mbtool.set_faceting_tol_tag(facet_tol.tolerance);

  for (const auto &p : j) {
    std::string step_file = p["filename"].get<std::string>();
    std::string material = p["material"].get<std::string>();

    std::cout << step_file << " : " << material << std::endl;
    std::vector<TopoDS_Shape> breps = step_to_breps(step_file);
    std::vector<std::string> emptyList;
    for (TopoDS_Shape shape : breps) {
      sew_and_facet2(shape, facet_tol, mbtool, emptyList, material);
    }
  }

  mbtool.write_geometry(h5m_file.c_str());
}
