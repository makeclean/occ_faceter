#ifndef STEP2BREPS_HH
#define STEP2BREPS_HH 1

#include <vector>
#include <string>
#include "TopoDS_Shape.hxx"

std::vector<TopoDS_Shape> step_to_breps(std::string step_file);

#endif // STEP2BREPS_HH
