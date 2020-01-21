#include "dagmc_faceter.hh"

#include <iostream>
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
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"

#include "TopLoc_Location.hxx"

#include "BRepOffsetAPI_Sewing.hxx"

#include "MBTool.hpp"

#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "NCollection_IndexedDataMap.hxx"

#include "BRepGProp.hxx"
#include "GProp_GProps.hxx"
#include "UniqueId/UniqueId.h" // this file also dep "UniqueId/half.hpp"
#include "read_metadata.hh"

typedef NCollection_IndexedDataMap<TopoDS_Face, moab::EntityHandle, TopTools_ShapeMapHasher> MapFaceToSurface;

struct TriangulationWithLocation {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation;
};

facet_data make_surface_facets(const TopoDS_Face &currentFace, const TriangulationWithLocation &facetData) {
  facet_data facets_for_moab;

  Handle(Poly_Triangulation) triangles = facetData.triangulation;
  const gp_Trsf &local_transform = facetData.loc;

  if (triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return facets_for_moab;
  } else {
    // retrieve facet data

    const TColgp_Array1OfPnt &nodes = triangles->Nodes();
    for (int i = nodes.Lower(); i <= nodes.Upper(); i++) {
      Standard_Real x, y, z;
      nodes(i).Coord(x, y, z);
      local_transform.Transforms(x, y, z);
      std::array<double, 3> coordinates;
      coordinates[0] = x, coordinates[1] = y, coordinates[2] = z;
      facets_for_moab.coords.push_back(coordinates);
    }
    // copy the facet_data
    std::array<int, 3> conn;
    const Poly_Array1OfTriangle &tris = triangles->Triangles();
    //     std::cout << "Face has " << tris.Length() << " triangles" << std::endl;
    for (int i = tris.Lower(); i <= tris.Upper(); i++) {
      // get the node indexes for this triangle
      const Poly_Triangle &tri = tris(i);
      tri.Get(conn[0], conn[1], conn[2]);

      facets_for_moab.connectivity.push_back(conn);
    }
    return facets_for_moab;
  }
}

// make the edge facets
edge_data make_edge_facets(const TopoDS_Edge &currentEdge,
                           const TriangulationWithLocation &facetData) {

  edge_data edges_for_moab;

  if (facetData.triangulation.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
  } else {
    // get the faceting for the edge
    Handle(Poly_PolygonOnTriangulation) edges =
        BRep_Tool::PolygonOnTriangulation(currentEdge, facetData.triangulation, facetData.loc);

    // convert TColStd_Array1OfInteger to std::vector<int>
    std::vector<int> &conn = edges_for_moab.connectivity;
    const TColStd_Array1OfInteger &lines = edges->Nodes();
    for (int i = lines.Lower(); i <= lines.Upper(); i++) {
      conn.push_back(lines(i));
    }
  }
  return edges_for_moab;
}

struct surface_data {
  facet_data facets;
  std::vector<edge_data> edge_collection;
};

surface_data get_facets_for_face(const TopoDS_Face &currentFace) {
  surface_data surface;

  // get the triangulation for the current face
  TriangulationWithLocation data;
  data.triangulation = BRep_Tool::Triangulation(currentFace, data.loc);

  // make facets for current face
  surface.facets = make_surface_facets(currentFace, data);

  TopTools_IndexedMapOfShape edges;
  TopExp::MapShapes(currentFace, TopAbs_EDGE, edges);
  for (int i = 1; i <= edges.Extent(); i++) {
    const TopoDS_Edge &currentEdge = TopoDS::Edge(edges(i));
    // make the edge facets
    edge_data edges = make_edge_facets(currentEdge, data);
    surface.edge_collection.push_back(edges);
  }
  return surface;
}

// Use BRepMesh_IncrementalMesh to make the triangulation
void perform_faceting(const TopoDS_Face &face, float facet_tol) {
  // This constructor calls Perform()
  BRepMesh_IncrementalMesh facets(face, facet_tol, false, 0.5);
}

std::uint64_t calculate_unique_id(const TopoDS_Shape &shape) {
  GProp_GProps v_props;
  BRepGProp::VolumeProperties(shape, v_props);
  return PPP::Utilities::geometryUniqueID(v_props.Mass(),
                                          {v_props.CentreOfMass().X(),
                                           v_props.CentreOfMass().Y(),
                                           v_props.CentreOfMass().Z()});
}

void facet_all_volumes(const TopTools_HSequenceOfShape &shape_list,
                       float facet_tol, MBTool &mbtool, MaterialsMap &mat_map) {
  int count = shape_list.Length();

  std::vector<TopoDS_Face> uniqueFaces;
  MapFaceToSurface surfaceMap;

  // list unique faces, create empty surfaces, and build surfaceMap
  for (int i = 1; i <= count; i++) {
    const TopoDS_Shape &shape = shape_list.Value(i);
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
      // Important note: For the surface map, face equivalence is defined
      // by TopoDS_Shape::IsSame(), which ignores the orientation.
      if (surfaceMap.Contains(face))
        continue;

      uniqueFaces.push_back(face);
      moab::EntityHandle surface;
      mbtool.make_new_surface(surface);
      surfaceMap.Add(face, surface);
    }
  }

  // do the hard work
  // (a range based for loop doesn't seem to work with OpenMP)
#pragma omp parallel for
  for (int i = 0; i < uniqueFaces.size(); i++) {
    perform_faceting(uniqueFaces[i], facet_tol);
  }

  // add facets (and edges) to surfaces
  for (MapFaceToSurface::Iterator it(surfaceMap); it.More(); it.Next()) {
    const TopoDS_Face &face = it.Key();
    moab::EntityHandle surface = it.Value();
    surface_data data = get_facets_for_face(face);
    mbtool.add_facets_and_curves_to_surface(surface, data.facets, data.edge_collection);
  }

  // build a list of volumes for each material
  std::map<std::string, std::vector<moab::EntityHandle>> material_volumes;

  // create volumes and add surfaces
  for (int i = 1; i <= count; i++) {
    const TopoDS_Shape &shape = shape_list.Value(i);

    moab::EntityHandle vol;
    mbtool.make_new_volume(vol);

    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
      moab::EntityHandle surface = surfaceMap.FindFromKey(face);
      int sense = face.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
      mbtool.add_surface_to_volume(surface, vol, sense);
    }

    // update map from material name to volumes
    if (!mat_map.empty()) {
      std::uint64_t uniqueID = calculate_unique_id(shape);
      std::string material = mat_map[uniqueID];
      if (!material.empty()) {
        material_volumes[material].push_back(vol);
      } else  {
        std::cout << "No material found for ID: " << uniqueID << std::endl;
      }
    }
  }

  // create material groups
  for (auto &pair : material_volumes) {
    std::string material = pair.first;
    std::vector<moab::EntityHandle> &volumes = pair.second;
    mbtool.add_group(material, volumes);
  }
}

void sew_shapes(const TopoDS_Shape &shape, TopTools_HSequenceOfShape &sewed_shapes) {
  if (shape.ShapeType() == TopAbs_COMPOUND) {
    // decend and get children
    for (TopoDS_Iterator it = TopoDS_Iterator(shape); it.More(); it.Next()) {
      sew_shapes(it.Value(), sewed_shapes);
    }
  } else if (shape.ShapeType() == TopAbs_SOLID) {
    // sew together all the curves
    BRepOffsetAPI_Sewing sew;
    sew.Add(shape);
    sew.Perform();

    // insert into the list
    sewed_shapes.Append(sew.SewedShape());
  } else {
    std::cout << "Unknown shape type " << shape.ShapeType() << std::endl;
  }
}

void sew_and_facet(TopoDS_Shape &shape, float facet_tol, MBTool &mbtool,
                   MaterialsMap &mat_map) {
  TopTools_HSequenceOfShape shape_list;
  sew_shapes(shape, shape_list);
  std::cout << "Instanciated " << shape_list.Length() << " items from file" << std::endl;

  facet_all_volumes(shape_list, facet_tol, mbtool, mat_map);
}

void dagmc_faceter(std::string brep_file, float facet_tol, std::string h5m_file) {
  TopoDS_Shape shape;
  BRep_Builder builder;
  BRepTools::Read(shape, brep_file.c_str(), builder);

  std::string json_file = brep_file.substr(0, brep_file.length() - 5) + "_metadata.json";

  MaterialsMap materials_map;
  read_metadata(json_file, materials_map);

  MBTool mbtool;
  mbtool.set_tags();
  sew_and_facet(shape, facet_tol, mbtool, materials_map);
  mbtool.write_geometry(h5m_file.c_str());
}
