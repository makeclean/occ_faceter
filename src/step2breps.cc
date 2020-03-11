#include "step2breps.hh"
#include <iostream>

#include "STEPControl_Reader.hxx"

std::vector<TopoDS_Shape> step_to_breps(std::string step_file) {
  STEPControl_Reader step;
  step.ReadFile(step_file.c_str());
  step.PrintCheckLoad(false, IFSelect_ListByItem);
  step.ClearShapes();
  int count = step.NbRootsForTransfer();
  std::cout << "STEP NbRootsForTransfer: " << count << std::endl;
  std::cout << "STEP TransferRoots: " << step.TransferRoots() << std::endl;
  std::cout << "STEP NbShapes: " << step.NbShapes() << std::endl;

  std::vector<TopoDS_Shape> breps;
  for (int i = 1; i <= count; i++) {
    bool ok = step.TransferRoot(i);
    step.PrintCheckTransfer(false, IFSelect_CountByItem);
    if (ok) {
      breps.push_back(step.Shape(i));
    }
  }
  return breps;
}
