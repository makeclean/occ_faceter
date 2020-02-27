#include "steps2h5m.hh"

#include <iostream>
#include <fstream>

#include "STEPControl_Reader.hxx"
#include "brep_faceter.hh"
#include "MBTool.hpp"
#include "step2breps.hh"
#include "UniqueId/json.hpp"

using json = nlohmann::json;

void steps2h5m(std::string input_file, double facet_tol, std::string h5m_file) {

  json j = json::parse(std::ifstream(input_file));

  MBTool mbtool;
  mbtool.set_tags();

  for (const auto &p : j) {
    std::string step_file = p["filename"].get<std::string>();
    std::string material = p["material"].get<std::string>();

    std::cout << step_file << " : " << material << std::endl;
    std::vector<TopoDS_Shape> breps = step_to_breps(step_file);
    MaterialsMap emptyMap;
    for (TopoDS_Shape shape : breps) {
      sew_and_facet(shape, facet_tol, mbtool, emptyMap, material);
    }
  }

  mbtool.write_geometry(h5m_file.c_str());
}
