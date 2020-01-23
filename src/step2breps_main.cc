#include <iostream>
#include "BRepTools.hxx"
#include "step2breps.hh"

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage: step2breps STEP_FILE OUTPUT_FILES_PREFIX" << std::endl;
    return 1;
  }

  std::string step_file(argv[1]);
  std::string brep_file_prefix(argv[2]);

  std::vector<TopoDS_Shape> breps = step_to_breps(step_file);
  for (int i = 0; i < breps.size(); i++) {
    std::string path = brep_file_prefix + std::to_string(i) + ".brep";
    BRepTools::Write(breps[i], path.c_str());
    std::cout << i << " ";
  }

  std::cout << std::endl;
  return 0;
}
