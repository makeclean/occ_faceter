#include <iostream>

#include "STEPControl_Reader.hxx"
#include "TopoDS_Shape.hxx"
#include "BRepTools.hxx"

int step2breps(std::string step_file, std::string brep_file_prefix) {
  STEPControl_Reader *step = new STEPControl_Reader();

  step->ReadFile(step_file.c_str());
  step->PrintCheckLoad(false, IFSelect_ListByItem);
  step->ClearShapes();
  int count = step->NbRootsForTransfer();
  int r_count = step->TransferRoots();
  std::cout << "count: " << count << std::endl;
  std::cout << "r_count: " << r_count << std::endl;
  std::cout << "n_shapes: " << step->NbShapes() << std::endl;

  for (int i = 1; i <= count; i++) {
    bool ok = step->TransferRoot(i);
    step->PrintCheckTransfer(false, IFSelect_CountByItem);
    if (ok) {
      TopoDS_Shape shape = step->Shape(i);
      std::string path = brep_file_prefix + std::to_string(i) + ".brep";
      BRepTools::Write(shape, path.c_str());
      std::cout << i << " ";
    }
  }

  std::cout << std::endl;

  delete step;
  return count;
}

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cout << "Usage: step2breps STEP_FILE OUTPUT_FILES_PREFIX" << std::endl;
    return 1;
  }

  std::string step_file(argv[1]);
  std::string brep_file_prefix(argv[2]);
  step2breps(step_file, brep_file_prefix);

  return 0;
}
