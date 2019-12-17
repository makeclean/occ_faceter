#include <iostream>
#include <array>
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

#include "XSControl_WorkSession.hxx"
#include "Interface_InterfaceModel.hxx"
#include "StepRepr_Representation.hxx"
#include "TCollection_HAsciiString.hxx"

#include "MBTool.hpp"

#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "NCollection_IndexedDataMap.hxx"

typedef NCollection_IndexedDataMap<TopoDS_Face, moab::EntityHandle, TopTools_ShapeMapHasher> MapFaceToSurface;

MBTool *mbtool = new MBTool();

float facet_tol = 0.;

struct FaceterData {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation;
};

// get the triangulation for the current face
FaceterData get_triangulation(TopoDS_Face currentFace) {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangles = BRep_Tool::Triangulation(currentFace, loc);
  FaceterData data;
  data.loc = loc;
  data.triangulation = triangles;
  return data;
}

facet_data make_surface_facets(TopoDS_Face currentFace, FaceterData facetData) {
  facet_data facets_for_moab;

  Handle(Poly_Triangulation) triangles = facetData.triangulation;
  TopLoc_Location loc = facetData.loc;

  if (triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return facets_for_moab;
  } else {
    // retrieve facet data

    gp_Trsf local_transform = loc;
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
      Poly_Triangle tri = tris(i);
      tri.Get(conn[0], conn[1], conn[2]);

      facets_for_moab.connectivity.push_back(conn);
    }
    return facets_for_moab;
  }
}

// make the edge facets
edge_data make_edge_facets(TopoDS_Face currentFace, TopoDS_Edge currentEdge, FaceterData facetData) {

  edge_data edges_for_moab;
  Handle(Poly_Triangulation) triangles = facetData.triangulation;
  TopLoc_Location loc = facetData.loc;

  // get the faceting for the edge
  Handle(Poly_PolygonOnTriangulation) edges = BRep_Tool::PolygonOnTriangulation(currentEdge, triangles, loc);

  if (triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return edges_for_moab;
  } else {
    // retrieve facet data
    std::vector<int> conn;
    const TColStd_Array1OfInteger &lines = edges->Nodes();
    for (int i = lines.Lower(); i <= lines.Upper(); i++) {
      conn.push_back(lines(i));
    }
    edges_for_moab.connectivity = conn;
    return edges_for_moab;
  }
}

struct surface_data {
  facet_data facets;
  std::vector<edge_data> edge_collection;
};

surface_data get_facets_for_face(TopoDS_Face currentFace) {
  surface_data surface;

  // get the triangulation for the current face
  FaceterData data = get_triangulation(currentFace);
  // make facets for current face
  surface.facets = make_surface_facets(currentFace, data);

  TopTools_IndexedMapOfShape edges;
  TopExp::MapShapes(currentFace, TopAbs_EDGE, edges);
  for (int i = 1; i <= edges.Extent(); i++) {
    TopoDS_Edge currentEdge = TopoDS::Edge(edges(i));
    // make the edge facets
    edge_data edges = make_edge_facets(currentFace, currentEdge, data);
    surface.edge_collection.push_back(edges);
  }
  return surface;
}

// Use BRepMesh_IncrementalMesh to make the triangulation
void perform_faceting(TopoDS_Face face) {
  // This constructor calls Perform()
  BRepMesh_IncrementalMesh facets(face, facet_tol, false, 0.5);
}

void facet_all_volumes(Handle_TopTools_HSequenceOfShape shape_list) {
  int count = shape_list->Length();

  std::vector<TopoDS_Face> uniqueFaces;
  MapFaceToSurface surfaceMap;

  // list unique faces, create empty surfaces, and build surfaceMap
  for (int i = 1; i <= count; i++) {
    TopoDS_Shape shape = shape_list->Value(i);
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      TopoDS_Face face = TopoDS::Face(ex.Current());
      // Important note: For the surface map, face equivalence is defined
      // by TopoDS_Shape::IsSame(), which ignores the orientation.
      if (surfaceMap.Contains(face))
        continue;

      uniqueFaces.push_back(face);
      moab::EntityHandle surface;
      mbtool->make_new_surface(surface);
      surfaceMap.Add(face, surface);
    }
  }

  // do the hard work
  // (a range based for loop doesn't seem to work with OpenMP)
#pragma omp parallel for
  for (int i = 0; i < uniqueFaces.size(); i++) {
    perform_faceting(uniqueFaces[i]);
  }

  // add facets (and edges) to surfaces
  for (MapFaceToSurface::Iterator it(surfaceMap); it.More(); it.Next()) {
    TopoDS_Face face = it.Key();
    moab::EntityHandle surface = it.Value();
    surface_data data = get_facets_for_face(face);
    mbtool->add_facets_and_curves_to_surface(surface, data.facets, data.edge_collection);
  }

  // create volumes and add surfaces
  for (int i = 1; i <= count; i++) {
    TopoDS_Shape shape = shape_list->Value(i);

    moab::EntityHandle vol;
    mbtool->make_new_volume(vol);

    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      TopoDS_Face face = TopoDS::Face(ex.Current());
      moab::EntityHandle surface = surfaceMap.FindFromKey(face);
      int sense = face.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
      mbtool->add_surface_to_volume(surface, vol, sense);
    }
  }
}

STEPControl_Reader *step;

//dump all labels
void dump_labels() {
  const Handle(XSControl_WorkSession) &theSession = step->WS();
  const Handle(Interface_InterfaceModel) &theModel = theSession->Model();

  Standard_Integer nb = theModel->NbEntities();
  for (Standard_Integer i = 1; i <= nb; i++) {

    Handle(StepRepr_Representation) ent =
        Handle(StepRepr_Representation)::DownCast(theModel->Value(i));
    if (ent.IsNull()) continue;
    if (ent->Name().IsNull()) continue;
    std::cout << ent->Name()->ToCString() << std::endl;
  }
  return;
}

void sew_and_append(TopoDS_Shape shape, Handle(TopTools_HSequenceOfShape) & shape_list) {
  // sew together all the curves
  BRepOffsetAPI_Sewing(1.0e-06, Standard_True);
  BRepOffsetAPI_Sewing sew;
  sew.Add(shape);
  sew.Perform();
  shape = sew.SewedShape();
  shape_list->Append(shape);
  return;
}

int main(int argc, char *argv[]) {

  std::string cad_file(argv[1]);
  std::string ftol(argv[2]);
  facet_tol = std::stod(ftol);
  std::string filename(argv[3]);

  moab::ErrorCode rval = mbtool->set_tags();

  step = new STEPControl_Reader();
  step->ReadFile(cad_file.c_str());

  step->PrintCheckLoad(false, IFSelect_ListByItem);
  step->ClearShapes();
  int count = step->NbRootsForTransfer();
  TopoDS_Shape shape;

  int r_count = step->TransferRoots();
  std::cout << "r_count: " << r_count << std::endl;
  std::cout << "n_shapes: " << step->NbShapes() << std::endl;

  Handle(TopTools_HSequenceOfShape) shape_list = new TopTools_HSequenceOfShape;
  std::cout << count << std::endl;
  for (int i = 1; i <= count; i++) {
    bool ok = step->TransferRoot(i);
    step->PrintCheckTransfer(false, IFSelect_CountByItem);
    if (ok) {
      shape = step->Shape(i);
      // if its a compound decend and get its children
      if (shape.ShapeType() == 0) {
        TopoDS_Iterator it = TopoDS_Iterator(shape);
        while (it.More()) {
          // if we are a volume
          if (it.Value().ShapeType() == 2) {
            sew_and_append(it.Value(), shape_list);
          } else {
            std::cout << "Unknown shape type " << it.Value().ShapeType() << std::endl;
          }
          it.Next();
        }
      } else if (shape.ShapeType() == 2) {
        // we are a normal volume insert into
        // the list
        sew_and_append(shape, shape_list);
      } else {
        std::cout << "Unknown shape" << std::endl;
      }
    } else {
      std::cout << "Couldnt read shape " << std::endl;
    }
  }
  std::cout << "Instanciated " << shape_list->Length() << " items from file" << std::endl;

  std::cout << count << " entities read from file " << std::endl;

  //dump_labels();

  facet_all_volumes(shape_list);
  mbtool->write_geometry(filename.c_str());
  delete mbtool;
  return 0;
}
