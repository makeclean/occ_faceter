#include <iostream>
#include <fstream>

#include "STEPControl_Reader.hxx"
#include "dagmc_faceter.hh"
#include "read_metadata.hh"
#include "MBTool.hpp"
#include "UniqueId/json.hpp"
#include "step2breps.hh"

using json = nlohmann::json;

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: steps2h5m JSON_FILE TOLERANCE H5M_FILE" << std::endl;
    return 1;
  }

  std::ifstream json_stream(argv[1]);
  float facet_tol = std::stof(argv[2]);
  std::string h5m_file(argv[3]);

  json j = json::parse(json_stream);

  MBTool mbtool;
  mbtool.set_tags();

  for (const auto &p : j) {
    std::string step_file = p["filename"].get<std::string>();
    std::string material = p["material"].get<std::string>();

    // add "mat:"" prefix to non-empty materials, unless it's already there
    if (!material.empty() && material.rfind("mat:", 0) != 0) {
      material = "mat:" + material;
    }

    std::cout << step_file << " : " << material << std::endl;
    std::vector<TopoDS_Shape> breps = step_to_breps(step_file);
    MaterialsMap emptyMap;
    for (TopoDS_Shape shape : breps) {
      sew_and_facet(shape, facet_tol, mbtool, emptyMap, material);
    }
  }

  mbtool.write_geometry(h5m_file.c_str());

  return 0;
}
