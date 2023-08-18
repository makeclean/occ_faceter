#include "brep_faceter.hh"

#include <iostream>
#include <algorithm>
#include <array>
#include <map>
#include <unordered_map>
#include <omp.h>

#include "BRep_Tool.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "TColStd_Array1OfInteger.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_PolygonOnTriangulation.hxx"

#include "Interface_Static.hxx"
#include "STEPControl_Reader.hxx"

#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Edge.hxx"

#include "TopTools_HSequenceOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"

#include "TopLoc_Location.hxx"

#include "BRepOffsetAPI_Sewing.hxx"

#include "MBTool.hpp"

#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "NCollection_IndexedDataMap.hxx"

// Terminology:
//
// Open Cascade Shape - corresponding moab meshset
//  [Material] - Group
//  Solid      - Volume
//  Face       - Surface
//  Edge       - Curve
//  Vertex     - Vertex
//
// Triangles from triangulation / faceted surfaces:
//  Triangles and Nodes (not "verticies")
//  MBEdge is two nodes, PolyEdge is multiple nodes
//
// "creating" a meshset includes adding it to the respective map
// "populating" a meshset can include creating and adding children and contents

typedef NCollection_IndexedDataMap<TopoDS_Face, moab::EntityHandle, TopTools_ShapeMapHasher> MapFaceToSurface;
typedef NCollection_IndexedDataMap<TopoDS_Edge, moab::EntityHandle, TopTools_ShapeMapHasher> MapEdgeToCurve;
typedef NCollection_IndexedDataMap<TopoDS_Vertex, moab::EntityHandle, TopTools_ShapeMapHasher> MapVertexToMeshset;

class BrepFaceter {
public:
  BrepFaceter(MBTool &mbt) : mbtool(mbt),
                             degenerate_triangle_count(0),
                             surface_without_facet_count(0) {}

  entity_vector facet(const TopTools_HSequenceOfShape &shape_list,
                      const FacetingTolerance &facet_tol);

private:
  MBTool &mbtool;
  MapFaceToSurface surfaceMap;
  MapEdgeToCurve edgeMap;
  MapVertexToMeshset vertexMap;
  entity_vector volumesList;

  int degenerate_triangle_count;
  int surface_without_facet_count;

  void create_surfaces(const TopTools_HSequenceOfShape &shape_list);
  void perform_faceting(const FacetingTolerance &facet_tol);
  void populate_all_surfaces();
  void create_volumes_and_add_children(const TopTools_HSequenceOfShape &shape_list);

  void create_surface_nodes(entity_vector &nodes, const Poly_Triangulation &triangulation, const TopLoc_Location &location);
  void create_surface_triangles(entity_vector &triangles, moab::EntityHandle surface, const Poly_Triangulation &triangulation, const entity_vector &nodes);

  void populate_vertex(moab::EntityHandle meshset, const TopoDS_Vertex &currentVertex) {
    // create and add contents
    double x, y, z;
    BRep_Tool::Pnt(currentVertex).Coord().Coord(x, y, z);
    moab::EntityHandle node = mbtool.find_or_create_node({x, y, z});

    mbtool.add_entity(meshset, node);
  }

  void populate_curve(moab::EntityHandle curve,
                      const TopoDS_Edge &currentEdge,
                      const Handle(Poly_Triangulation) &triangulation,
                      const TopLoc_Location &location,
                      const entity_vector &surfaceNodes) {
    // get the faceting for the edge
    Handle(Poly_PolygonOnTriangulation) edges =
        BRep_Tool::PolygonOnTriangulation(currentEdge, triangulation, location);

    if (edges.IsNull()) {
      std::cerr << "Warning: Unexpected null edges." << std::endl;
      return;
    }

    const TColStd_Array1OfInteger &lines = edges->Nodes();
    if (lines.Length() < 2) {
      std::cerr << "Warning: Attempting to build empty curve." << std::endl;
      return;
    }

    entity_vector mbedge_entities;
    entity_vector node_entities;

    auto occ_edge_it = lines.cbegin();

    // subtract one because OCC uses one based indexing
    moab::EntityHandle prev = surfaceNodes.at(*occ_edge_it - 1);
    node_entities.push_back(prev);
    occ_edge_it++;

    for (; occ_edge_it != lines.cend(); occ_edge_it++) {
      moab::EntityHandle node = surfaceNodes.at(*occ_edge_it - 1);
      moab::EntityHandle edge = mbtool.create_edge({prev, node});

      node_entities.push_back(node);
      mbedge_entities.push_back(edge);
      prev = node;
    }

    // if curve is closed, remove duplicate vertex
    if (node_entities.front() == node_entities.back()) {
      node_entities.pop_back();
    }

    // add nodes and mbedges to curve
    mbtool.add_entities(curve, node_entities);
    mbtool.add_entities(curve, mbedge_entities);

    // create and populate children
    for (TopExp_Explorer explorer(currentEdge, TopAbs_VERTEX); explorer.More(); explorer.Next()) {
      const TopoDS_Vertex &currentVertex = TopoDS::Vertex(explorer.Current());

      moab::EntityHandle meshset;
      if (!vertexMap.FindFromKey(currentVertex, meshset)) {
        meshset = mbtool.make_new_vertex();
        vertexMap.Add(currentVertex, meshset);

        populate_vertex(meshset, currentVertex);
      }

      mbtool.add_child_to_parent(meshset, curve);
    }
  }

  void populate_surface(moab::EntityHandle surface, const TopoDS_Face &face) {
    // get the triangulation for the current face
    TopLoc_Location location;
    Handle(Poly_Triangulation) triangulation = BRep_Tool::Triangulation(face, location);

    if (triangulation.IsNull() || triangulation->NbNodes() < 1) {
      surface_without_facet_count += 1;
      return;
    }

    // add contents to surface
    entity_vector nodes;
    create_surface_nodes(nodes, triangulation, location);
    mbtool.add_entities(surface, nodes);

    entity_vector triangles;
    create_surface_triangles(triangles, surface, triangulation, nodes);
    mbtool.add_entities(surface, triangles);

    // recursively create, populate and add childen
    for (TopExp_Explorer edges(face, TopAbs_EDGE); edges.More(); edges.Next()) {
      const TopoDS_Edge &currentEdge = TopoDS::Edge(edges.Current());
      moab::EntityHandle curve;
      if (!edgeMap.FindFromKey(currentEdge, curve)) {
        curve = mbtool.make_new_curve();
        edgeMap.Add(currentEdge, curve);
        populate_curve(curve, currentEdge, triangulation, location, nodes);
      }
      int sense = currentEdge.Orientation() != face.Orientation() ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
      mbtool.add_child_to_parent(curve, surface, sense);
    }
  }

  void populate_volume(moab::EntityHandle vol, const TopoDS_Shape &shape) {
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
      moab::EntityHandle surface = surfaceMap.FindFromKey(face);
      int sense = face.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
      mbtool.add_child_to_parent(surface, vol, sense);
    }
  }
};

void BrepFaceter::create_surfaces(const TopTools_HSequenceOfShape &shape_list) {
  // list unique faces, create empty surfaces, and build surface map

  // Important note: For the maps, edge/face equivalence is defined
  // by TopoDS_Shape::IsSame(), which ignores the orientation.
  for (const TopoDS_Shape &shape : shape_list) {
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
      if (surfaceMap.Contains(face))
        continue;

      surfaceMap.Add(face, mbtool.make_new_surface());
    }
  }
}

void BrepFaceter::create_surface_nodes(entity_vector &nodes,
                                       const Poly_Triangulation &triangulation,
                                       const TopLoc_Location &location) {
  const gp_Trsf &local_transform = location;
  // retrieve facet data
  for (int i = 1; i <= triangulation.NbNodes(); i++) {
    Standard_Real x, y, z;
    triangulation.Node(i).Coord(x, y, z);
    local_transform.Transforms(x, y, z);
    nodes.push_back(mbtool.find_or_create_node({x, y, z}));
  }
}

void BrepFaceter::create_surface_triangles(entity_vector &triangles,
                                           moab::EntityHandle surface,
                                           const Poly_Triangulation &triangulation,
                                           const entity_vector &nodes) {
  //     std::cout << "Face has " << tris.Length() << " triangles" << std::endl;
  for (int i = 1; i <= triangulation.NbTriangles(); i++) {
    // get the node indexes for this triangle
    const Poly_Triangle &tri = triangulation.Triangle(i);

    int a, b, c;
    tri.Get(a, b, c);
    // subtract one because OCC uses one based indexing
    std::array<moab::EntityHandle, 3> connections = {
        nodes.at(a - 1),
        nodes.at(b - 1),
        nodes.at(c - 1),
    };
    if (connections[2] == connections[1] ||
        connections[1] == connections[0] ||
        connections[2] == connections[0]) {
      degenerate_triangle_count += 1;
    } else {
      triangles.push_back(mbtool.create_triangle(connections));
    }
  }
}

void BrepFaceter::perform_faceting(const FacetingTolerance &facet_tol) {

#pragma omp parallel for
  for (int i = 1; i <= surfaceMap.Extent(); i++) {
    // This constructor calls Perform() to mutate the face adding triangulation
    // that can be used by the following serial code
    BRepMesh_IncrementalMesh(surfaceMap.FindKey(i), facet_tol.tolerance, facet_tol.is_relative, 0.5);
  }
}

void BrepFaceter::populate_all_surfaces() {
  // Note: surface meshsets have actually been created early, but they are
  // populated here, and edge and vertex meshsets are created here.

  for (MapFaceToSurface::Iterator it(surfaceMap); it.More(); it.Next()) {
    const TopoDS_Face &face = it.Key();
    moab::EntityHandle surface = it.Value();

    populate_surface(surface, face);
  }

  // Instead of outputting "No facets for surface" or "degenerate triangle
  // not created" every time it occurs, we count occurences and then output
  // a single warning message saying how many time the issue has occured
  // (following an issue with many lines of unhelpful output when processing
  // one example geometry).
  if (surface_without_facet_count > 0) {
    std::cerr << "Warning: " << surface_without_facet_count
              << " surfaces found without facets." << std::endl;
  }

  if (degenerate_triangle_count > 0) {
    std::cerr << "Warning: " << degenerate_triangle_count
              << " degenerate triangles have been ignored." << std::endl;
  }
}

void BrepFaceter::create_volumes_and_add_children(const TopTools_HSequenceOfShape &shape_list) {
  for (const TopoDS_Shape &shape : shape_list) {
    moab::EntityHandle vol = mbtool.make_new_volume();
    volumesList.push_back(vol);

    populate_volume(vol, shape);
  }
}

entity_vector BrepFaceter::facet(const TopTools_HSequenceOfShape &shape_list,
                                 const FacetingTolerance &facet_tol) {
  // build surfaceMap (with early creations of surface meshsets) so we have a
  // set of unique faces for faceting
  create_surfaces(shape_list);

  perform_faceting(facet_tol);
  populate_all_surfaces();

  // Note: Not much further "population" of volumes required, since they have
  // no contents (they do have children) and their children have already been
  // created and populated.
  create_volumes_and_add_children(shape_list);

  return volumesList;
}

static std::string add_mat_prefix(const std::string &material) {
  // add "mat:"" prefix, unless it's already there
  if (!material.empty() && material.rfind("mat:", 0) != 0) {
    return "mat:" + material;
  }
  return material;
}

void add_materials(MBTool &mbtool, const entity_vector &volumes,
                   const std::vector<std::string> &mat_list) {
  // build a list of volumes for each material
  std::map<std::string, entity_vector> material_volumes;

  if (mat_list.size() < volumes.size()) {
    std::cerr << "Error: Material list is too short." << std::endl;
    return;
  }

  // keep track of where we are in the material list
  auto next_material = mat_list.begin();

  // update map from material name to volumes
  for (moab::EntityHandle vol : volumes) {
    const std::string &material = *(next_material++);
    material_volumes[material].push_back(vol);
  }

  // create material groups
  for (auto &pair : material_volumes) {
    mbtool.add_group(add_mat_prefix(pair.first), pair.second);
  }
}

void add_single_material(MBTool &mbtool, const entity_vector &volumes,
                         const std::string &single_material) {
  mbtool.add_group(add_mat_prefix(single_material), volumes);
}

entity_vector sew_and_facet2(TopoDS_Shape &shape, const FacetingTolerance &facet_tol, MBTool &mbtool) {
  TopTools_HSequenceOfShape shape_list;
  for (TopExp_Explorer solids(shape, TopAbs_SOLID); solids.More(); solids.Next()) {
    // sew together all the curves
    BRepOffsetAPI_Sewing sew;
    sew.Add(solids.Current());
    sew.Perform();

    // insert into the list
    shape_list.Append(sew.SewedShape());
  }
  std::cout << "Instanciated " << shape_list.Length() << " items from file" << std::endl;

  BrepFaceter bf(mbtool);
  return bf.facet(shape_list, facet_tol);
}

void read_materials_list(std::string text_file, std::vector<std::string> &mat_list) {
  std::ifstream text_stream(text_file);
  if (text_stream.fail()) {
    std::cerr << "Warning: Failed to read file " << text_file << std::endl;
  } else {
    std::string line;
    while (std::getline(text_stream, line)) {
      mat_list.push_back(line);
    }
  }
}

void brep_faceter(std::string brep_file, std::string materials_list_file,
                  const FacetingTolerance &facet_tol, std::string h5m_file,
                  bool add_mat_ids, double scale_factor) {
  TopoDS_Shape shape;
  BRep_Builder builder;
  BRepTools::Read(shape, brep_file.c_str(), builder);

  std::vector<std::string> materials_list;
  read_materials_list(materials_list_file, materials_list);

  MBTool mbtool;
  // TODO: review use of GEOMETRY_RESABS
  mbtool.set_faceting_tol_tag(facet_tol.tolerance);
  mbtool.set_scale_factor(scale_factor);
  entity_vector volumes = sew_and_facet2(shape, facet_tol, mbtool);
  add_materials(mbtool, volumes, materials_list);

  if (add_mat_ids)
    mbtool.add_mat_ids();

  mbtool.gather_ents();
  mbtool.write_geometry(h5m_file.c_str());
}
