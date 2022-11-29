// moab includes
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"

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


#include "vertex_inserter.hh"

// system includes
#include <fstream>
#include <iostream>

// cgal typedefs
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef std::vector<K::Point_3> Polyline_type;
typedef std::list< Polyline_type > Polylines;
typedef CGAL::AABB_halfedge_graph_segment_primitive<Mesh> HGSP;
typedef CGAL::AABB_traits<K, HGSP>    AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>  AABB_tree;
typedef K::Point_3 Vertex;

//typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::Polyhedron_3<K,CGAL::Polyhedron_items_with_id_3>  Polyhedron;
//typedef CGAL::Polyhedron_items_with_id_3<K> Polyhedron; 
typedef Polyhedron::HalfedgeDS HalfedgeDS;
typedef CGAL::Nef_polyhedron_3<K> Nef_Polyhedron;
//typedef Polyhedron::Vertex Vertex;
//typedef Polyhedron::Vertex_iterator Vertex_iterator;
//typedef Polyhedron::Facet_iterator Facet_iterator;

moab::Core *MBI = NULL;

namespace moab {
class GeomTopoTool;
} 

int surfID = 0;
int volID = 0;

const char geom_categories[][CATEGORY_TAG_SIZE] = {"Vertex\0",
						   "Curve\0",
						   "Surface\0",
						   "Volume\0",
						   "Group\0"};

  moab::Tag geometry_dimension_tag, id_tag;
  moab::Tag faceting_tol_tag, geometry_resabs_tag;
  moab::Tag category_tag;
  moab::Tag vol_id_tag, surf_id_tag; // tags for triangles for plotting
  moab::Tag name_tag;
  moab::Tag mat_id_tag;

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

  void generateSurfaceSegmentation(Polyhedron p) {
    // create a property-map for SDF values
    typedef std::map<Polyhedron::Facet_const_handle, double> Facet_double_map;
    Facet_double_map internal_sdf_map;
    boost::associative_property_map<Facet_double_map> sdf_property_map(internal_sdf_map);
    // compute SDF values using default parameters for number of rays, and cone angle
    CGAL::sdf_values(p, sdf_property_map);
    // create a property-map for segment-ids
    typedef std::map<Polyhedron::Facet_const_handle, std::size_t> Facet_int_map;
    Facet_int_map internal_segment_map;
    boost::associative_property_map<Facet_int_map> segment_property_map(internal_segment_map);
    // segment the mesh using default parameters for number of levels, and smoothing lambda
    // Any other scalar values can be used instead of using SDF values computed using the CGAL function
    std::size_t number_of_segments = CGAL::segmentation_from_sdf_values(p, sdf_property_map, segment_property_map);
    std::cout << "Number of segments: " << number_of_segments << std::endl;
    // print segment-ids
    for(Polyhedron::Facet_const_iterator facet_it = p.facets_begin();
        facet_it != p.facets_end(); ++facet_it) {
        // ids are between [0, number_of_segments -1]
        std::cout << segment_property_map[facet_it] << " ";
    }
    std::cout << std::endl;

    const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
    const double smoothing_lambda = 0.3;  // importance of surface features, suggested to be in-between [0,1]
    // Note that we can use the same SDF values (sdf_property_map) over and over again for segmentation.
    // This feature is relevant for segmenting the mesh several times with different parameters.
    CGAL::segmentation_from_sdf_values(
      p, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda);
     for(Polyhedron::Facet_const_iterator facet_it = p.facets_begin();
        facet_it != p.facets_end(); ++facet_it) {
        // ids are between [0, number_of_segments -1]
        std::cout << segment_property_map[facet_it] << " ";
     }
     std::cout << std::endl;

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
      std::cout << *it << std::endl;
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
    /*

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
    */
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

  // get the CGAL version of geometry
  void makeCGALGeometryPoly() {
   std::map<int, moab::EntityHandle>::iterator it;
    for ( it = volmap.begin() ; it != volmap.end() ; ++it ) {
      Polyhedron p = getCGALMeshAsPolyhedron(it->first);
      geometry_poly[it->first] = p;
    }   
  }

  // subtract vol1 from vol2
  void subtract_volumes_poly(int vol1, int vol2){
    Polyhedron p1 = geometry_poly[vol1];
    Polyhedron p2 = geometry_poly[vol2];

    moab::ErrorCode rval = moab::MB_FAILURE;

    Nef_Polyhedron np1(p1);
    Nef_Polyhedron np2(p2);
    
    Nef_Polyhedron up2 = np1 - np2;
    Polyhedron update2;
    up2.convert_to_polyhedron(update2);
    up2.clear();

    rval = nvm->addPolyhedron(update2,1);
    rval = nvm->addPolyhedron(p2,2);

    nvm->write("new_geom.h5m");

    return;
  }
  
  /*
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
  */

  /*
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
  */
  
  private:
  NewVolsMOAB *nvm;
  moab::Tag geometry_dimension_tag, id_tag;
  std::map<int, moab::EntityHandle> volmap;
  std::map<int, Mesh> geometry;
  std::map<int, Polyhedron> geometry_poly;
};

moab::ErrorCode add_surface_to_volume(moab::EntityHandle surface,
          moab::EntityHandle volume, int sense) {
  moab::ErrorCode rval;
  rval = MBI->add_parent_child(volume, surface);
  MB_CHK_ERR(rval);
  moab::GeomTopoTool *geom_tool = new moab::GeomTopoTool(MBI);
  rval = geom_tool->set_sense(surface, volume, sense);
  delete geom_tool;
  return rval;
}

//  makes a new surface in moab
moab::ErrorCode make_new_surface(moab::EntityHandle &surface) {
  surfID++;
  //  moab::EntityHandle surface;
  // std::cout << "Created new surface " << surfID << std::endl;
  
  moab::ErrorCode rval = MBI->create_meshset(moab::MESHSET_ORDERED, surface);
  // set the id tag
  rval = MBI->tag_set_data(id_tag,&surface,1,&surfID);
  // set the dim tag
  int dim = 2;
  rval = MBI->tag_set_data(geometry_dimension_tag,&surface,1,&dim);
  rval = MBI->tag_set_data(category_tag,&surface,1,&geom_categories[dim]);
  return moab::MB_SUCCESS;  
}

// make a new volume meshset 
moab::ErrorCode make_new_volume(moab::EntityHandle &volume) {
  volID++;

  // std::cout << "Created new volume " << volID << std::endl;
  
  // make a new volume set
  moab::ErrorCode rval = MBI->create_meshset(moab::MESHSET_ORDERED,volume);
  // set the id tag
  rval = MBI->tag_set_data(id_tag,&volume,1,&volID);
  // set the dim tag
  int dim = 3;
  rval = MBI->tag_set_data(geometry_dimension_tag,&volume,1,&dim);
  rval = MBI->tag_set_data(category_tag,&volume,1,&geom_categories[dim]);
  return rval;
}

// 
int main(int argc, char* argv[]) {
  MBI = new moab::Core();
    moab::ErrorCode rval = moab::MB_FAILURE;

  rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Couldnt get geom dim tag");
    rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Couldnt get id tag");

    rval = MBI->tag_get_handle("VOL_ID", 1, moab::MB_TYPE_INTEGER, vol_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Couldnt get vol_id tag");

    rval = MBI->tag_get_handle("SURF_ID", 1, moab::MB_TYPE_INTEGER, surf_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Couldnt get surf_id tag");
    
    rval = MBI->tag_get_handle("FACETING_TOL", 1, moab::MB_TYPE_DOUBLE, faceting_tol_tag,
			       moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Error creating faceting_tol_tag");
    rval = MBI->tag_get_handle("GEOMETRY_RESABS", 1, moab::MB_TYPE_DOUBLE, 
			       geometry_resabs_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Error creating geometry_resabs_tag");
    
    rval = MBI->tag_get_handle(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, moab::MB_TYPE_OPAQUE, 
    			       category_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Error creating category_tag");

    rval = MBI->tag_get_handle(NAME_TAG_NAME, NAME_TAG_SIZE, moab::MB_TYPE_OPAQUE,
                               name_tag, moab::MB_TAG_SPARSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Error creating name_tag");

    rval = MBI->tag_get_handle("MatID", 1, moab::MB_TYPE_INTEGER,
                               mat_id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
    MB_CHK_SET_ERR(rval, "Error creating mat_id_tag");

  moab::EntityHandle input_set;
  rval = MBI->create_meshset(moab::MESHSET_SET, input_set); 
  MB_CHK_SET_ERR(rval, "failed to create meshset");
  std::cout << "Loading input file..." << std::endl;
  rval = MBI->create_meshset(moab::MESHSET_SET, input_set); 
  MB_CHK_SET_ERR(rval, "failed to create meshset");

  moab::EntityHandle s1,v1, s2,v2;
  rval = make_new_surface(s1);
  rval = make_new_volume(v1);
  rval = add_surface_to_volume(s1,v1,1);
  rval = make_new_surface(s2);
  rval = make_new_volume(v2);
  rval = add_surface_to_volume(s2,v2,1);

  rval = MBI->load_file("../../vol1.stl", &s1);

  rval = MBI->load_file("../../vol2.stl", &s2);
  
  std::cout << rval << std::endl;

  rval = MBI->write_file("intermediate.h5m");

  // make the cgal geometry
  MOABInterface *cgal = new MOABInterface();
  cgal->makeCGALGeometryPoly();
  cgal->subtract_volumes_poly(1,2);
  return 0;
}
