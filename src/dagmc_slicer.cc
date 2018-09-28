#include "dagmc_slicer.hh"

// system includes
#include <fstream>
#include <iostream>

MOABInterface::MOABInterface(moab::Core *moab){
  MBI = moab;
  setupInterface();
  getVolumeIDList();
}
  
MOABInterface::~MOABInterface(){
}
  
// setup 
void MOABInterface::setupInterface() {
  moab::ErrorCode rval;
  rval = MBI->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, moab::MB_TYPE_INTEGER, geometry_dimension_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Couldnt get geom dim tag");
  rval = MBI->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, moab::MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE | moab::MB_TAG_CREAT);
  MB_CHK_SET_ERR_RET(rval, "Couldnt get id tag");
}

// get a list of all the volumes in the problem
void MOABInterface::getVolumeIDList() {
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
moab::Range MOABInterface::getChildTriangles(int vol_id) {
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
Mesh MOABInterface::makeCGALMesh(moab::Range triangles) {
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
Mesh MOABInterface::getCGALMesh(int vol_id) {
  moab::ErrorCode rval;
  moab::Range triangles = getChildTriangles(vol_id);
  Mesh mesh = makeCGALMesh(triangles);
  return mesh;
}

// slice through the whole geometry
std::map<int,Polylines> MOABInterface::sliceGeometry(double dir[3], double offset) {
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
void MOABInterface::makeCGALGeometry() {
  std::map<int, moab::EntityHandle>::iterator it;
  for ( it = volmap.begin() ; it != volmap.end() ; ++it ) {
    Mesh mesh = getCGALMesh(it->first);
    geometry[it->first] = mesh;
  }
}

