#define BOOST_PARAMETER_MAX_ARITY 12

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
//#include "CGAL/Polygon_mesh_processing/corefinement.h"
#include "CGAL/Point_3.h"
#include "CGAL/predicates_on_points_3.h"
#include "CGAL/Polyhedron_3.h"
#include "CGAL/Polyhedron_incremental_builder_3.h"
#include "CGAL/Nef_polyhedron_3.h"
#include <CGAL/Polyhedron_items_with_id_3.h> 
// mesh segmentation
#include <CGAL/mesh_segmentation.h>
#include <CGAL/property_map.h>
// 3d meshing
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

#include "vertex_inserter.hh"

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

//typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedron_3<K,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
//typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;

//typedef CGAL::Polyhedron_items_with_id_3<K> Polyhedron; 
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef CGAL::Nef_polyhedron_3<K> Nef_Polyhedron;
//typedef Polyhedron::Vertex Vertex;
//typedef Polyhedron::Vertex_iterator Vertex_iterator;
//typedef Polyhedron::Facet_iterator Facet_iterator;

moab::Core *MBI = NULL;

using namespace CGAL::parameters;

const char geom_categories[][CATEGORY_TAG_SIZE] = {"Vertex\0",
						   "Curve\0",
						   "Surface\0",
						   "Volume\0",
						   "Group\0"};

// class to store the new volume information
class NewVolsMOAB {
  public:
  // constructor
  NewVolsMOAB(){
    MOAB = new moab::Core();
    vi = new VertexInserter::VertexInserter(MOAB);
    setup_tags();
  }

  void setup_tags() {
    moab::ErrorCode rval = moab::MB_FAILURE;
    rval = MOAB->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");

    rval = MOAB->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");
    
    rval = MOAB->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
		  	     moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating faceting_tol_tag");

    rval = MOAB->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
		  	     geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating geometry_resabs_tag");
    
    rval = MOAB->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
	  		     category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR_RET(rval, "Error creating category_tag");

    rval = MOAB->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
	  		     name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
  
    return;
  }

  // set the surface tags
  void set_surface_tags(int id, moab::EntityHandle surface_set) {
    const void* id_val = &id;
    moab::ErrorCode rval;

    int val_s = 2;
    const void* val_surf = &val_s;
    rval = MOAB->tag_set_data(category_tag,&surface_set,1,&geom_categories[2]);
    rval = MOAB->tag_set_data(geometry_dimension_tag,&surface_set,1,val_surf);
    rval = MOAB->tag_set_data(id_tag,&surface_set,1,id_val);
    return;
  }

  // set the volume tags
  void set_volume_tags(int id, moab::EntityHandle volume_set) {
    const void* id_val = &id;
    moab::ErrorCode rval;

    int val_v = 3;
    const void* val_vol = &val_v;
    rval = MOAB->tag_set_data(category_tag,&volume_set,1,&geom_categories[3]);
    rval = MOAB->tag_set_data(geometry_dimension_tag,&volume_set,1,val_vol);
    rval = MOAB->tag_set_data(id_tag,&volume_set,1,id_val);
    return;
  }

  // get trangles and verts from polyhedron
  moab::ErrorCode getTrisFromPoly(Polyhedron p, moab::Range &triangles) {
    std::cout << p.size_of_facets() << std::endl;
    moab::ErrorCode rval;
    if(p.size_of_facets() == 0 ) return moab::MB_FAILURE;
    //
    std::unordered_map<int,moab::EntityHandle> index_handle_map;
    int idx = 0;
    for ( Polyhedron::Vertex_iterator v = p.vertices_begin() ; 
          v != p.vertices_end() ; ++v ) {
        v->id() = idx++;
        CGAL::Point_3<K> point = v->point();
        double x = CGAL::to_double(point.x());
        double y = CGAL::to_double(point.y());
        double z = CGAL::to_double(point.z());
        std::array<double,3> coords = {x,y,z};
        moab::EntityHandle handle;
        rval = vi->insert_vertex(coords,handle);
        index_handle_map[v->id()] = handle;
    } 

    idx = 0;
    // loop over the facets in the shape
    for ( Polyhedron::Facet_iterator  f = p.facets_begin() ; 
          f != p.facets_end() ; ++f ) {
      f->id() = idx++;
      Polyhedron::Halfedge_around_facet_circulator circ = f->facet_begin(); 
      std::vector<moab::EntityHandle> connectivity;
      do { 
        int id = circ->vertex()->id(); 
        moab::EntityHandle handle = index_handle_map[id];
        connectivity.push_back(handle);
      } while ( ++circ != f->facet_begin());
      // make a new triangle
      assert(connectivity.size() == 3);
      moab::EntityHandle triangle;
      rval = MOAB->create_element(moab::MBTRI,&connectivity[0],3,triangle); 
      triangles.insert(triangle);
    }

    return moab::MB_SUCCESS;
  }

  // add polyedron with id id
  moab::ErrorCode addPolyhedron(Polyhedron p, int id) {
    moab::ErrorCode rval = moab::MB_FAILURE;
    moab::Range surface;

    // todo - figure out why we get 0's
    // generateSurfaceSegmentation(p);

    rval = getTrisFromPoly(p,surface);

    // tag the surface set
    moab::EntityHandle surface_set;
    rval = MOAB->create_meshset(moab::MESHSET_SET,surface_set);
    set_surface_tags(id,surface_set);
 
    // tag the volume set
    moab::EntityHandle volume_set;
    rval = MOAB->create_meshset(moab::MESHSET_SET,volume_set);
    set_volume_tags(id,volume_set);

    // add parent child links
    rval = MOAB->add_parent_child(volume_set,surface_set);
    // add the triangles to the surface set 
    rval = MOAB->add_entities(surface_set,surface); 
  
  }

  void write(std::string filename) {
    moab::ErrorCode rval = moab::MB_FAILURE;
    rval = MOAB->write_mesh(filename.c_str());
    return;
  }

  // destructor
  ~NewVolsMOAB(){
    delete MOAB;
  }
  moab::Core *MOAB;
  VertexInserter::VertexInserter *vi;
  // tags
  moab::Tag geometry_dimension_tag;
  moab::Tag id_tag;
  moab::Tag faceting_tol_tag;
  moab::Tag geometry_resabs_tag;
  moab::Tag category_tag;
  moab::Tag name_tag;
};


class MOABInterface {
  public:
  MOABInterface(){
    nvm = new NewVolsMOAB();
    setupInterface();
    getVolumeIDList();
  }
  ~MOABInterface(){
    delete nvm;
  }

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class Build_triangle : public CGAL::Modifier_base<HDS> {  
  public:
  Build_triangle(moab::Core *mbi,
		 const moab::Range verts,
		 const moab::Range tris) {
    MBI = mbi;
    vertices = verts;
    triangles = tris;
  }
  
    void operator()( HDS& hds) {
      // Postcondition: hds is a valid polyhedral surface.

      CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
      B.begin_surface( vertices.size(), triangles.size(),
		       triangles.size()*3, 0 );
      
      typedef typename HDS::Vertex Vertex;
      typedef typename Vertex::Point Point;

      // add all the vertices
      std::map<moab::EntityHandle,int> cgal_verts;
      double pos[3];
      int idx = 0;
      moab::ErrorCode rval = moab::MB_FAILURE;
      for ( moab::EntityHandle vertex : vertices) {
	      rval = MBI->get_coords(&(vertex),1,&pos[0]);
	      B.add_vertex(Point(pos[0],pos[1],pos[2]));
	      cgal_verts[vertex] = idx;
	      idx++;
      }

      // loop over the triangles and add them
      //      B.begin_surface();
      std::vector<moab::EntityHandle> verts;
      for ( moab::EntityHandle triangle : triangles ) {
	      B.begin_facet();
	      rval = MBI->get_connectivity(&(triangle), 1, verts);
	      B.add_vertex_to_facet( cgal_verts[verts[0]]);
	      B.add_vertex_to_facet( cgal_verts[verts[1]]);
	      B.add_vertex_to_facet( cgal_verts[verts[2]]);
        B.end_facet();
      }
      B.end_surface();
    }
  
  private: 
    moab::Range triangles;
    moab::Range vertices;
    moab::Core *MBI;
};
  
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

      std::vector<Mesh::Vertex_index> points;
      points.push_back(point1);
      points.push_back(point2);
      points.push_back(point3);
      
      mesh.add_face(points);
    }
    return mesh;
  }

  // make a CGAL polyhedron object from a bag of triangles
  Polyhedron makeCGALPolyhedron(moab::Range triangles) {
    moab::Range vertices;
    moab::ErrorCode rval = moab::MB_FAILURE;
    rval = MBI->get_vertices(triangles, vertices);
    // now get make an indexable list
    moab::Range::iterator it;

    Polyhedron p;
    Build_triangle<HalfedgeDS> triangle(MBI,vertices,triangles);
    p.delegate(triangle);
    return p;
  }

  // for a given dagmc volume get a mesh
  Mesh getCGALMesh(int vol_id) {
    moab::ErrorCode rval;
    moab::Range triangles = getChildTriangles(vol_id);
    Mesh mesh = makeCGALMesh(triangles);
    return mesh;
  }

  // for a given dagmc volume get a mesh
  Polyhedron getCGALMeshAsPolyhedron(int vol_id) {
    moab::ErrorCode rval;
    moab::Range triangles = getChildTriangles(vol_id);
    Polyhedron p = makeCGALPolyhedron(triangles);
    return p;
  }

  // make a cgal version of the dagmc geoemtry
  void makeCGALGeometry() {
    std::map<int, moab::EntityHandle>::iterator it;
    for ( it = volmap.begin() ; it != volmap.end() ; ++it ) {
      Mesh mesh = getCGALMesh(it->first);
      geometry[it->first] = mesh;
    }
  }

  // get the CGAL version of geometry
  void makeCGALGeometryPoly() {
   std::map<int, moab::EntityHandle>::iterator it;
    for ( it = volmap.begin() ; it != volmap.end() ; ++it ) {
      Polyhedron p = getCGALMeshAsPolyhedron(it->first);
      geometry_poly[it->first] = p;
    }   
  }

  Polyhedron getPolygon(int id) { return geometry_poly[id];}
  
  private:
  NewVolsMOAB *nvm;
  moab::Tag geometry_dimension_tag, id_tag;
  std::map<int, moab::EntityHandle> volmap;
  std::map<int, Mesh> geometry;
  std::map<int, Polyhedron> geometry_poly;
};

class CGALMesher {

  #ifdef CGAL_CONCURRENT_MESH_3
    typedef CGAL::Parallel_tag Concurrency_tag;
  #else
    typedef CGAL::Sequential_tag Concurrency_tag;
  #endif
  
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;
  //typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain; 
  //typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
  // Triangulation
  // Triangulation
  typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
  // Criteria
  //typedef CGAL::Mesh_triangulation_3<Mesh_domain,CGAL::Default,Concurrency_tag>::type Tr;
  //typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_index> C3t3;
  // Criteria
  typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;

  typedef typename C3t3::Cells_in_complex_iterator Cell_iterator;

  typedef typename Tr::Vertex_handle Vertex_handle;

  // To avoid verbose function and named parameters call
  public:
  CGALMesher() {
    MOAB = new moab::Core();
    VI = new VertexInserter::VertexInserter(MOAB);
  }
  ~CGALMesher();

  void generate_mesh(const Polyhedron p) {

    Mesh_domain domain(p);

    //domain.detect_features();
    
    // Mesh criteria (no cell_size set)
    //Mesh_criteria criteria(facet_angle=25, facet_size=0.15, facet_distance=0.008,
    //                       cell_radius_edge_ratio=3);
    Mesh_criteria criteria(cell_radius_edge_ratio=3,cell_size=0.3);
    
    // Mesh generation
    C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
					no_perturb(),
					no_exude() );
    
    // Mesh criteria (no cell_size set)
    moab::ErrorCode rval;
    
    // loop over 
    std::cout << c3t3.number_of_cells_in_complex() << std::endl;
    for( Cell_iterator cit = c3t3.cells_in_complex_begin();
         cit != c3t3.cells_in_complex_end();  ++cit) {
      // loop over the vertices
      std::vector<moab::EntityHandle>  entities(4); // 4 verts for a tet

      for ( int i = 0 ; i < 4 ; i++ ) {
	const Vertex_handle& vh = cit->vertex(i);
	double x = CGAL::to_double(vh->point().x());
	double y = CGAL::to_double(vh->point().y());
	double z = CGAL::to_double(vh->point().z());                        
	// make the vertex
	std::array<double,3> coords {x,y,z};
	rval = VI->insert_vertex(coords,entities[i]);
      }
      moab::EntityHandle tet;
      rval = MOAB->create_element(moab::MBTET,&(entities[0]),4,tet);
      // make the tet
    }
    rval = MOAB->write_mesh("3d_mesh.h5m");
    return;
  }
  
  moab::Core *MOAB;
  VertexInserter::VertexInserter *VI;
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
  cgal->makeCGALGeometryPoly();
  CGALMesher *mshr = new CGALMesher();
  mshr->generate_mesh(cgal->getPolygon(1));
  return 0;
}
