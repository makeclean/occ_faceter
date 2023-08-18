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

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#ifndef DAGMC_SLICER_HH 
#define DAGMC_SLICER_HH 1

// cgal typedefs
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef std::vector<K::Point_3> Polyline_type;
typedef std::list< Polyline_type > Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
typedef K::Point_3 Vertex;

class MOABInterface {
  public:
    MOABInterface(moab::Core *moab);
   ~MOABInterface();

  void setupInterface();
  void getVolumeIDList();
  moab::Range getChildTriangles(int vol_id);
  Mesh getCGALMesh(int vol_id);
  std::vector<Polylines> sliceGeometry(double dir[3], double offset);
  std::map<int,Polylines> sliceGeometryByID(double dir[3], double offset);
  void makeCGALGeometry();
  Mesh makeCGALMesh(moab::Range triangles);
  
  private:
  moab::Core *MBI;
  moab::Tag geometry_dimension_tag, id_tag;
  std::map<int, moab::EntityHandle> volmap;
  std::map<int, Mesh> geometry;  
};

#endif
