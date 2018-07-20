// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

// cgal includes
#include "CGAL/Exact_predicates_exact_constructions_kernel.h"
#include "CGAL/Surface_mesh.h"
#include "CGAL/AABB_halfedge_graph_segment_primitive.h"
#include "CGAL/AABB_tree.h"
#include "CGAL/AABB_traits.h"
#include "CGAL/Polygon_mesh_slicer.h"
#include "CGAL/Polygon_mesh_processing/corefinement.h"
#include "CGAL/Point_3.h"
#include "CGAL/predicates_on_points_3.h"

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
      /*
      switch (CGAL::orientation(point1,point2,point3)) {
        case CGAL::LEFTTURN: std::cout << "Left Turn\n"; break;
        case CGAL::RIGHTURN: std::cout << "Right Turn\n"; break;
        case CGAL::COLLINEAR: std::cout << "No Turn\n"; break;	
      }
      */
      //mesh.add_face(point1,point2,point3);
      std::vector<Mesh::Vertex_index> points;
      points.push_back(point1);
      points.push_back(point2);
      points.push_back(point3);
      
      mesh.add_face(points);
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

  // given two volumes - subtract them
  void subtract_volumes() {
    std::cout << "Subtract volumes" << std::endl;
    Mesh mesh1 = geometry[1];
    Mesh mesh2 = geometry[2];    

    // check that the meshes are filled
    if(!CGAL::is_triangle_mesh(mesh1)) { std::cout << "Mesh not filled" << std::endl; exit(0);}
    if(!CGAL::is_triangle_mesh(mesh2)) { std::cout << "Mesh not filled" << std::endl; exit(0);}
    
    Mesh out;
    //bool valid = true;
    std::ofstream output1("refined1.off"); output1 << mesh1 ; output1.close();
    std::ofstream output2("refined2.off"); output2 << mesh2 ; output2.close();

    std::cout << "corefining" << std::endl;
    CGAL::Polygon_mesh_processing::corefine(mesh1,mesh2);
    //    std::ofstream output1("refined1.off"); output1 << mesh1 ; output1.close();
    //std::ofstream output2("refined2.off"); output2 << mesh2 ; output2.close();
    std::cout << "unioning" << std::endl;
    bool valid = CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh1,mesh2,out);
    std::cout << "Valid: " << valid << std::endl;
    if(!valid) {
      std::cout << "failed to compute intersection" << std::endl;
    } else {
      std::ofstream output("intersection.off");
      output << out;
    }
    return;
  }

  // given two volumes - subtract them
  void subtract_volumes_mesh() {
    std::cout << "Subtract volumes" << std::endl;
    Mesh mesh1;// = geometry[1];
    Mesh mesh2;// = geometry[2];

    std::ifstream m1("refined1.off");
    std::ifstream m2("refined2.off");

    //std::cout << !m1 << " " << !m2 << std::endl;
    if ( !m1 ) { std::cout << "file doesnt exist" << std::endl;}
    if ( !m2 ) { std::cout << "file doesnt exist" << std::endl;}
    if ( !(m1 >> mesh1)) {
      std::cerr << "Not a valid input file." << std::endl;
    }

    if ( !(m2 >> mesh2) ) {
      std::cerr << "Not a valid input file." << std::endl;
    }
    
    
    m1.close();
    m2.close();

    //CGAL::triangulate_polyhedron(mesh1);
    //CGAL::triangulate_polyhedron(mesh2);
    
    
    std::cout << CGAL::is_triangle_mesh(mesh1) << " ";
    std::cout << CGAL::is_triangle_mesh(mesh2) << std::endl;
        
    Mesh out;
    //bool valid = true;
    CGAL::Polygon_mesh_processing::corefine(mesh1,mesh2);
    std::ofstream output1("refined_1.off"); output1 << mesh1 ; output1.close();
    std::ofstream output2("refined_2.off"); output2 << mesh2 ; output2.close();
    std::cout << "woof" << std::endl;
    return;
    bool valid = CGAL::Polygon_mesh_processing::corefine_and_compute_union(mesh1,mesh2,out);
    std::cout << "Valid: " << valid << std::endl;
    if(!valid) {
      std::cout << "failed to compute intersection" << std::endl;
    } else {
      std::ofstream output("intersection.off");
      output << out;
    }
    return;
  }
  
  private:
  moab::Tag geometry_dimension_tag, id_tag;
  std::map<int, moab::EntityHandle> volmap;
  std::map<int, Mesh> geometry;
};



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
  cgal->subtract_volumes_mesh();
  
  return 0;
}
