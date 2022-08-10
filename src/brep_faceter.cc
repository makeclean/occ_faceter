#include "brep_faceter.hh"

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
typedef NCollection_IndexedDataMap<TopoDS_Edge, moab::EntityHandle, TopTools_ShapeMapHasher> MapEdgeToCurve;

struct TriangulationWithLocation {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation;
};

facet_data make_surface_facets(const TopoDS_Face &currentFace,
                               const TriangulationWithLocation &facetData) {
  facet_data facets_for_moab;

  Handle(Poly_Triangulation) triangles = facetData.triangulation;
  const gp_Trsf &local_transform = facetData.loc;

  // retrieve facet data
  for (int i = 1; i <= triangles->NbNodes(); i++) {
    Standard_Real x, y, z;
    triangles->Node(i).Coord(x, y, z);
    local_transform.Transforms(x, y, z);
    std::array<double, 3> coordinates;
    coordinates[0] = x, coordinates[1] = y, coordinates[2] = z;
    facets_for_moab.coords.push_back(coordinates);
  }
  // copy the facet_data
  std::array<int, 3> conn;
  //     std::cout << "Face has " << tris.Length() << " triangles" << std::endl;
  for (int i = 1; i <= triangles->NbTriangles(); i++) {
    // get the node indexes for this triangle
    const Poly_Triangle &tri = triangles->Triangle(i);
    tri.Get(conn[0], conn[1], conn[2]);

    facets_for_moab.connectivity.push_back(conn);
  }
  return facets_for_moab;
}

// make the edge facets
edge_data make_edge_facets(const TopoDS_Edge &currentEdge,
                           const TriangulationWithLocation &facetData) {
  edge_data edges_for_moab;

  if (facetData.triangulation.IsNull())
    return edges_for_moab;

  // get the faceting for the edge
  Handle(Poly_PolygonOnTriangulation) edges =
      BRep_Tool::PolygonOnTriangulation(currentEdge, facetData.triangulation, facetData.loc);

  if (edges.IsNull()) {
    std::cout << "Warning: Unexpected null edges." << std::endl;
    return edges_for_moab;
  }

  // convert TColStd_Array1OfInteger to std::vector<int>
  std::vector<int> &conn = edges_for_moab.connectivity;
  const TColStd_Array1OfInteger &lines = edges->Nodes();
  // TODO: Confirm this is correct when currentEdge.Orientation() == TopAbs_REVERSED
  //   - in /Cubit-plugin/.../DAGMCExportCommand.cpp the equivalent list of
  //     edge nodes is reversed for these edges.
  for (int i = lines.Lower(); i <= lines.Upper(); i++) {
    conn.push_back(lines(i));
  }
  return edges_for_moab;
}

// Use BRepMesh_IncrementalMesh to make the triangulation
void perform_faceting(const TopoDS_Face &face, const FacetingTolerance& facet_tol) {
  // This constructor calls Perform()
  BRepMesh_IncrementalMesh facets(face, facet_tol.tolerance, facet_tol.is_relative, 0.5);
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
                       const FacetingTolerance& facet_tol,
                       MBTool &mbtool, MaterialsMap &mat_map,
                       std::string single_material) {
  int count = shape_list.Length();

  std::vector<TopoDS_Face> uniqueFaces;
  MapFaceToSurface surfaceMap;
  MapEdgeToCurve edgeMap;

  // list unique faces, create empty surfaces, and build surface map

  // Important note: For the maps, edge/face equivalence is defined
  // by TopoDS_Shape::IsSame(), which ignores the orientation.
  for (int i = 1; i <= count; i++) {
    const TopoDS_Shape &shape = shape_list.Value(i);
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
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
  int n_surfaces_without_facets = 0;
  for (MapFaceToSurface::Iterator it(surfaceMap); it.More(); it.Next()) {
    const TopoDS_Face &face = it.Key();
    moab::EntityHandle surface = it.Value();

    // get the triangulation for the current face
    TriangulationWithLocation data;
    data.triangulation = BRep_Tool::Triangulation(face, data.loc);

    if (data.triangulation.IsNull() || data.triangulation->NbNodes() < 1) {
      n_surfaces_without_facets++;
    } else {
      // make facets for current face
      facet_data facets = make_surface_facets(face, data);
      facet_vertex_map f_vertex_map;
      moab::ErrorCode ret = mbtool.generate_facet_vertex_map(f_vertex_map, facets.coords);
      assert(ret == moab::MB_SUCCESS);
      mbtool.add_facets_to_surface(surface, facets.connectivity, f_vertex_map);

      // add curves to surface
      TopTools_IndexedMapOfShape edges;
      TopExp::MapShapes(face, TopAbs_EDGE, edges);
      for (int i = 1; i <= edges.Extent(); i++) {
        const TopoDS_Edge &currentEdge = TopoDS::Edge(edges(i));

        moab::EntityHandle curve;
        if (!edgeMap.FindFromKey(currentEdge, curve)) {
          mbtool.make_new_curve(curve);
          edgeMap.Add(currentEdge, curve);

          edge_data edges = make_edge_facets(currentEdge, data);
          mbtool.build_curve(curve, edges, f_vertex_map);
        }
        int sense = currentEdge.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
        mbtool.add_child_to_parent(curve, surface, sense);
      }
    }
  }

  if (n_surfaces_without_facets > 0) {
    std::cout << "Warning: " << n_surfaces_without_facets
      << " surfaces found without facets." << std::endl;
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
      mbtool.add_child_to_parent(surface, vol, sense);
    }

    // update map from material name to volumes

    // if single_material is set, then ignore mat_map
    if (!single_material.empty()) {
      material_volumes[single_material].push_back(vol);
    } else if (!mat_map.empty()) {
      std::uint64_t uniqueID = calculate_unique_id(shape);
      std::string material = mat_map[uniqueID];
      if (!material.empty()) {
        material_volumes[material].push_back(vol);
      } else {
        std::cout << "No material found for ID: " << uniqueID << std::endl;
      }
    }
  }

  // create material groups
  for (auto &pair : material_volumes) {
    std::string material = pair.first;
    std::vector<moab::EntityHandle> &volumes = pair.second;

    // add "mat:"" prefix, unless it's already there
    if (!material.empty() && material.rfind("mat:", 0) != 0) {
      material = "mat:" + material;
    }

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

void sew_and_facet(TopoDS_Shape &shape, const FacetingTolerance& facet_tol, MBTool &mbtool,
                   MaterialsMap &mat_map, std::string single_material) {
  TopTools_HSequenceOfShape shape_list;
  sew_shapes(shape, shape_list);
  std::cout << "Instanciated " << shape_list.Length() << " items from file" << std::endl;

  facet_all_volumes(shape_list, facet_tol, mbtool, mat_map, single_material);
}

void brep_faceter(std::string brep_file, std::string json_file,
                  const FacetingTolerance& facet_tol, std::string h5m_file,
                  bool add_mat_ids) {
  TopoDS_Shape shape;
  BRep_Builder builder;
  BRepTools::Read(shape, brep_file.c_str(), builder);

  MaterialsMap materials_map;
  read_metadata(json_file, materials_map);

  MBTool mbtool;
  // TODO: review use of GEOMETRY_RESABS
  mbtool.set_faceting_tol_tag(facet_tol.tolerance);
  sew_and_facet(shape, facet_tol, mbtool, materials_map);

  if (add_mat_ids)
    mbtool.add_mat_ids();

  mbtool.gather_ents();
  mbtool.write_geometry(h5m_file.c_str());
}
