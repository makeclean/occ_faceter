#include <iostream>
#include "brep_faceter.hh"
#include "steps2h5m.hh"
#include "moab/ProgOptions.hpp"

static bool has_ending(std::string s, std::string ending) {
  return s.rfind(ending) == s.length() - ending.length();
}

int main(int argc, char *argv[]) {
  ProgOptions po("occ_faceter: faceting a geometry and saving to a MOAB h5m file");

  std::string input_file;
  std::string output_file = "dagmc_not_watertight.h5m";
  double tolerance = 0.001;
  bool tol_is_absolute = false;
  bool add_mat_ids = false;

  po.addRequiredArg<std::string>("input_file",
                                 "Path to brep and json output from PPP, or json list of step files and materials",
                                 &input_file);

  po.addOpt<double>("tolerance,t", "Faceting tolerance (default " + std::to_string(tolerance) + ")", &tolerance);
  po.addOpt<void>("absolute_tol,a", "Treat faceting tolerance as absolute (default is relative)", &tol_is_absolute);
  po.addOpt<std::string>("output_file,o", "Path to output file (default "+ output_file + ")", &output_file);
  po.addOpt<void>("add_mat_ids,m", "Add MatID to every triangle (default " + std::string(add_mat_ids ? "true)" : "false)"), &add_mat_ids);

  po.parseCommandLine(argc, argv);

  if (!has_ending(output_file, ".h5m")) {
    std::cerr << "Warning: output file path should end with .h5m" << std::endl;
  }

  FacetingTolerance facet_tol(tolerance, tol_is_absolute);

  if (has_ending(input_file, ".json")) {
    steps2h5m(input_file, facet_tol, output_file);
  } else if (has_ending(input_file, ".brep")) {
    // expecting a json file with similar path to the brep file,
    // but with ".brep" replaced by "_metadata.json"
    std::string json_file = input_file.substr(0, input_file.length() - 5) + "_metadata.json";

    brep_faceter(input_file, json_file, facet_tol, output_file, add_mat_ids);
  } else {
    std::cerr << "Error: Path to input file must end with .json or .brep" << std::endl;
    return 1;
  }

  return 0;
}
