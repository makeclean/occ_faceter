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
  rval = MBI->create_meshset(moab::MESHSET_SET, input_set); 
  MB_CHK_SET_ERR(rval, "failed to create meshset");
  std::cout << "Loading input file..." << std::endl;
  rval = MBI->load_file("test.h5m", &input_set);  
  
  // make the cgal geometry
  MOABInterface *cgal = new MOABInterface();
  cgal->makeCGALGeometry(); 

  std::map<int,Polylines> slices;

  double dir[3] = {0.,1.,0.};
  slices = cgal->sliceGeometry(dir,0.);
  cgal->writeSlice(slices, "slice1.txt");
  dir[1] = 0.;
  dir[2] = 1.0;
  slices = cgal->sliceGeometry(dir,0.);
  cgal->writeSlice(slices, "slice2.txt");
  dir[2] = 0.0;
  dir[0] = 1.0;
  slices = cgal->sliceGeometry(dir,0.);
  cgal->writeSlice(slices, "slice3.txt");

  /*
  // find a given volume - or for all the triangles
  // rather than get all triangles at once, maybe need to go volume by
  // volume
  moab::Range triangles;
  rval = MBI->get_entities_by_type(0,moab::MBTRI, triangles);
  MB_CHK_SET_ERR(rval, "failed to find any triangles");
  // now get a list of vertices
  moab::Range vertices;
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
    //Vertex point = Vertex(pos[0],pos[1],pos[2]);
    //mesh.add_vertex(point);
    cgal_verts[*it] = idx;
  }

  // now for each triangle make the triangle
  for ( it = triangles.begin() ; it != triangles.end() ; ++it) {
    std::vector<moab::EntityHandle> verts;
    rval = MBI->get_connectivity(&(*it),1, verts);
    Mesh::Vertex_index point1 = cgal_verts[verts[0]];
    Mesh::Vertex_index point2 = cgal_verts[verts[1]];
    Mesh::Vertex_index point3 = cgal_verts[verts[2]];
			     
    // Mesh.add_face(v1,v2,v3); - add all the triangles;
    mesh.add_face(point1,point2,point3);     
  }
  
  CGAL::Polygon_mesh_slicer<Mesh, K> slicer(mesh);
  Polylines slice;
  std::cout << "slice 1 ..." << std::endl;
  slicer(K::Plane_3(0,0,1,0.0), std::back_inserter(slice));
  write_slice(slice,"slice1.txt");
  slice.clear();
  std::cout << "slice 2 ..." << std::endl;
  slicer(K::Plane_3(0,0,1,0.0), std::back_inserter(slice));
  write_slice(slice,"slice2.txt");
  slice.clear();
  std::cout << "slice 3 ..." << std::endl;
  slicer(K::Plane_3(0,1,0,0.0), std::back_inserter(slice));
  write_slice(slice,"slice3.txt");
  slice.clear();

  // write mesh to off
  std::ofstream outfile("mesh.off");
  outfile << mesh;
  */  
  return 0;

}
