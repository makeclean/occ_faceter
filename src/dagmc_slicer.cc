// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

// cgal includes
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Surface_mesh.h"
#include "CGAL/AABB_halfedge_graph_segment_primitive.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/Polygon_mesh_slicer.h"

// boost includes
#include "boost/program_options.hpp"

// system includes
#include <fstream>
#include <iostream>

// cgal typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef std::vector<K::Point_3> Polyline_type;
typedef std::list< Polyline_type > Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
typedef K::Point_3 Vertex;

moab::Core *MBI = NULL;

class MOABInterface {
  public:
  MOABInterface(){
    setupInterface();
    getVolumeIDList();
  }
  ~MOABInterface(){
  }
  
  // setup 
  void setupInterface() {
    moab::ErrorCode rval;
    rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");
    rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");
  }

  // get a list of all the volumes in the problem
  void getVolumeIDList() {
    int dim = 3;
    const void *val = &dim;
    // get the volumes
    moab::Range volumes;
    moab::ErrorCode rval = MBI->get_entities_by_type_and_tag
      (0,moab::MBENTITYSET ,&geometry_dimension_tag, &val, 1, volumes);
    moab::Range::iterator it;
    int id = 0;
    for ( it = volumes.begin() ; it != volumes.end() ; ++it ) {
      rval = MBI->tag_get_data(id_tag, &(*it), 1, &id);
      volmap[id] = *it;
    }
  }

  // get the triangles that belong to a given volume
  moab::Range getChildTriangles(int vol_id) {
    moab::EntityHandle vol = volmap[vol_id];
    moab::Range child_surfs;
    moab::ErrorCode rval;
    rval = MBI->get_child_meshsets(vol,child_surfs);
    moab::Range::iterator it;
    // for each child surface get triangles
    moab::Range triangle_set, triangles;
    // we should probably check for the merge tag here, and invert normals
    // if sense needs to be reversed
    for ( it = child_surfs.begin() ; it != child_surfs.end() ; ++it ) {
      rval = MBI->get_entities_by_type(*it, moab::MBTRI, triangles);
      triangle_set.merge(triangles);
    }
    return triangle_set;
  }

  // make a CGAL mesh object from a bag of triangles
  Mesh makeCGALMesh(moab::Range triangles) {
    moab::Range vertices;
    moab::ErrorCode rval = moab::MB_FAILURE;
    rval = MBI->get_vertices(triangles, vertices);
    // now get make an indexable list
    moab::Range::iterator it;
    
    std::map<moab::EntityHandle,Mesh::Vertex_index> cgal_verts;

    Mesh mesh;
    double pos[3];
    Mesh::Vertex_index idx;
    for ( it = vertices.begin() ; it != vertices.end() ; ++it ) {
      rval = MBI->get_coords(&(*it),1,&pos[0]);
      idx = mesh.add_vertex(Vertex(pos[0],pos[1],pos[2]));
      cgal_verts[*it] = idx;
    }

    // now for each triangle make the triangle
    for ( it = triangles.begin() ; it != triangles.end() ; ++it) {
      std::vector<moab::EntityHandle> verts;
      rval = MBI->get_connectivity(&(*it),1, verts);
      Mesh::Vertex_index point1 = cgal_verts[verts[0]];
      Mesh::Vertex_index point2 = cgal_verts[verts[1]];
      Mesh::Vertex_index point3 = cgal_verts[verts[2]];			     
      mesh.add_face(point1,point2,point3);     
    }
    return mesh;
  }
  
  // for a given dagmc volume get a mesh
  Mesh getCGALMesh(int vol_id) {
    moab::ErrorCode rval;
    moab::Range triangles = getChildTriangles(vol_id);
    Mesh mesh = makeCGALMesh(triangles);
    return mesh;
  }

  // slice through the whole geometry
  std::map<int,Polylines> sliceGeometry(double dir[3], double offset) {
    std::map<int, Mesh>::iterator it;
    std::map<int,Polylines> slices;
    for ( it = geometry.begin() ; it != geometry.end() ; ++it ) {
      CGAL::Polygon_mesh_slicer<Mesh, K> slicer(it->second);
      Polylines slice;
      slicer(K::Plane_3(dir[0],dir[1],dir[2],offset), std::back_inserter(slice));
      slices[it->first] = slice;
    }
    return slices;
  }

  // make a cgal version of the dagmc geoemtry
  void makeCGALGeometry() {
    std::map<int, moab::EntityHandle>::iterator it;
    for ( it = volmap.begin() ; it != volmap.end() ; ++it ) {
      Mesh mesh = getCGALMesh(it->first);
      geometry[it->first] = mesh;
    }
  }

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
  
  private:
  moab::Tag geometry_dimension_tag, id_tag;
  std::map<int, moab::EntityHandle> volmap;
  std::map<int, Mesh> geometry;
};


// write the slice
void write_slice(Polylines lines, std::string filename) {
  std::ofstream outfile(filename);
  Polylines::iterator it;
  Polyline_type::iterator jt;
  for ( it = lines.begin() ; it != lines.end() ; ++it ) {
    int idx = std::distance(lines.begin(), it);
    idx++;
    //    outfile << "set object " << idx << " polygon from \\"  << std::endl;
    for ( jt = it->begin() ; jt != it->end() ; ++jt) {
      //outfile << (*jt).x() << "," << (*jt).y() /* <<  "," << (*jt).z() */ << " to \\" << std::endl;
      outfile << (*jt).x() << " " << (*jt).y() /* <<  "," << (*jt).z() */ <<  std::endl;
    }
    outfile << (*(it->begin())).x() << " " << (*(it->begin())).y() /* << "," << (*(it->begin())).z() */ << std::endl;
    outfile << " " << std::endl;
  }
  outfile.close();
}

// 
int main(int argc, char* argv[]) {
  MBI = new moab::Core();
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
  MOABInterface *cgal = new MOABInterface();
  cgal->makeCGALGeometry(); 

  std::map<int,Polylines> slices;

  double dir[3] = {0.,1.,0.};
  slices = cgal->sliceGeometry(dir,slice_position[1]);
  cgal->writeSlice(slices, "slice1.txt");
  dir[1] = 0.;
  dir[2] = 1.0;
  slices = cgal->sliceGeometry(dir,slice_position[0]);
  cgal->writeSlice(slices, "slice2.txt");
  dir[2] = 0.0;
  dir[0] = 1.0;
  slices = cgal->sliceGeometry(dir,slice_position[2]);
  cgal->writeSlice(slices, "slice3.txt");

  return 0;

}
