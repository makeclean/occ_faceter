#include <iostream>
#include "dagmc_faceter.hh"

int main(int argc, char *argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: occ_faceter BREP_FILE TOLERANCE H5M_FILE" << std::endl;
    return 1;
  }

  std::string brep_file(argv[1]);
  float facet_tol = std::stof(argv[2]);
  std::string h5m_file(argv[3]);

  dagmc_faceter(brep_file, facet_tol, h5m_file);
  return 0;
}
