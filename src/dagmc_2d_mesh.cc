// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "vertex_inserter.hh"
#include "dagmc_slicer.hh"

#include <fstream>
#include "boost/program_options.hpp"

#include "CGAL/Polygon_mesh_slicer.h"
#include "CGAL/Constrained_Delaunay_triangulation_2.h"
#include "CGAL/Triangulation_face_base_with_info_2.h"
#include "CGAL/Polygon_2.h"

struct FaceInfo2 {
  FaceInfo2(){}
  int nesting_level;
  bool in_domain(){
    return nesting_level%2 == 1;
  }
};

typedef CGAL::Triangulation_vertex_base_2<K>                      Vb;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2,K>    Fbb;
typedef CGAL::Constrained_triangulation_face_base_2<K,Fbb>        Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>               TDS;
typedef CGAL::Exact_predicates_tag                                Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag>  CDT;
typedef CDT::Point Point;
typedef CGAL::Polygon_2<K> Polygon_2;
typedef CDT::Face_handle Face_handle;
typedef CDT::Vertex_handle Vertex_handle;

// stride through the slices and make a triangulation
void make_2d_triangulation(const std::map<int,Polylines> slices) {
  // loop over the slices
  std::vector<Polygon_2> polygons;

  moab::Core *MOAB = new moab::Core();
  VertexInserter::VertexInserter *vi = new VertexInserter::VertexInserter(MOAB);

  moab::ErrorCode rval;
  
  for ( std::pair<int,Polylines> boundary : slices ) {
    std::cout << boundary.first << std::endl;
    std::cout << boundary.second.size() << std::endl;

    // for each closed loop in the boundary 
    for ( Polyline_type segments : boundary.second ) {
      // for each point in a given loop
      Polygon_2 polygon;
      for ( K::Point_3 point : segments ) {
        polygon.push_back(Point(point.x(),point.y()));
      }
      std::cout << polygon << std::endl;
      // append the polygon to the vector
      polygons.push_back(polygon);
    }
  }

  CDT cdt;
  for (int i = 0 ; i < polygons.size() ; i++ ) {
    cdt.insert_constraint(polygons[i].vertices_begin(),
			  polygons[i].vertices_end(),
			  true);
  }

  CDT::Finite_faces_iterator t = cdt.faces_begin();

  int count = 0;
  for ( t = cdt.finite_faces_begin() ;
        t != cdt.finite_faces_end() ;
        t++ ) {

    count++;
    moab::EntityHandle triangle;
    std::array<moab::EntityHandle,3> connectivity;
    for ( int i = 0 ; i < 3 ; i++ ) {

      double x = CGAL::to_double(t->vertex(i)->point().hx());
      double y = CGAL::to_double(t->vertex(i)->point().hy());
        
      std::array<double,3> coords = {x,y,0.0};
      // make new moab vertex
      rval = vi->insert_vertex(coords,connectivity[i]);
    }
    // make a new triangle
    rval = MOAB->create_element(moab::MBTRI,&(connectivity[0]),3,triangle);
  }
  rval = MOAB->write_mesh("mesh.h5m");
  delete vi;
  delete MOAB;
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

  double dir[3] = {0.,0.,1.};
  slices = cgal->sliceGeometry(dir,slice_position[1]);

  make_2d_triangulation(slices);
  
  return 0;

}
