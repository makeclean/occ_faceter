// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "dagmc_slicer.hh"

#include <fstream>
#include "boost/program_options.hpp"

#include "CGAL/Polygon_mesh_slicer.h"

//
void writeSlice(std::map<int,Polylines> slices, std::string filename) {
    std::map<int,Polylines>::iterator slice_it;
    std::ofstream outfile(filename);      
    for ( slice_it = slices.begin() ; slice_it != slices.end() ; ++slice_it) {
      Polylines::iterator it;
      Polyline_type::iterator jt;
      for ( it = slice_it->second.begin() ; it != slice_it->second.end() ; ++it) {
	for ( jt = it->begin() ; jt != it->end() ; ++jt ) {
	  outfile << (*jt).x() << " " <<  (*jt).y() << " " << (*jt).z() << std::endl;
	}
	outfile << std::endl;
      }
    }
  }


int main(int argc, char* argv[]) {

  moab::Core *MBI = new moab::Core();
  
  moab::ErrorCode rval = moab::MB_FAILURE;
  moab::EntityHandle input_set;

  boost::program_options::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("input-file", boost::program_options::value<std::string>(), "input file")
    ("slice-position", boost::program_options::value< std::vector<double> >()->multitoken(), "slice_position")
    ;

  boost::program_options::positional_options_description p;
  p.add("input_file", -1);

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc,argv).
	    options(desc).positional(p).run(),vm);
  boost::program_options::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return -1;
  }

  std::string filename = "";
  if (vm.count("input-file")) {
    filename = vm["input-file"].as<std::string>();
  } else {
    std::cout << "input-file not specified" << std::endl;
    return -1;
  }

  std::vector<double> slice_position;
  if (vm.count("slice-position")) {
    slice_position = vm["slice-position"].as< std::vector<double> >();
    if(slice_position.size() != 3) {
      std::cout << "There needs to three elements" << std::endl;
      return -1;
    }
  } else {
    std::cout << "slice-position not specified" << std::endl;
    return -1;
  }

  std::cout << "open file " << filename << std::endl;
  std::cout << "slicing at " << slice_position[0] << " " << slice_position[1];
  std::cout << " " << slice_position[2] << std::endl;
  
  rval = MBI->create_meshset(moab::MESHSET_SET, input_set); 
  MB_CHK_SET_ERR(rval, "failed to create meshset");
  std::cout << "Loading input file..." << std::endl;
  rval = MBI->load_file(filename.c_str(), &input_set);  
  
  // make the cgal geometry
  MOABInterface *cgal = new MOABInterface(MBI);
  cgal->makeCGALGeometry(); 

  std::map<int,Polylines> slices;

  double dir[3] = {0.,1.,0.};
  slices = cgal->sliceGeometry(dir,slice_position[1]);
  writeSlice(slices, "slice1.txt");
  dir[1] = 0.;
  dir[2] = 1.0;
  slices = cgal->sliceGeometry(dir,slice_position[0]);
  writeSlice(slices, "slice2.txt");
  dir[2] = 0.0;
  dir[0] = 1.0;
  slices = cgal->sliceGeometry(dir,slice_position[2]);
  writeSlice(slices, "slice3.txt");

  return 0;

}
