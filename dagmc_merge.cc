#include <iostream>
#include "dagmc_topology.hpp"

// 
int main (int argc, char* argv[]) {
  // create topology check class
  DAGMCTopology *DT = new DAGMCTopology();
  // load the file
  moab::ErrorCode rval = moab::MB_FAILURE;
  rval = DT->load_file(std::string(argv[1]));
  // perform_merge
  rval = DT->perform_merge();
  // save the file
  rval = DT->save_file(std::string(argv[2]));
  
  return 0;
}
