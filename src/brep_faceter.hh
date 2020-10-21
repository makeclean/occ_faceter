#ifndef BREP_FACETER_HH
#define BREP_FACETER_HH 1

#include <iostream>
#include <array>
#include <map>
#include <unordered_map>
#include <omp.h>
#include <string>

#include "GProp_GProps.hxx"
#include "UniqueId/UniqueId.h" // this file also dep "UniqueId/half.hpp"
#include "read_metadata.hh"

#include "BRepOffsetAPI_Sewing.hxx"
#include "BRep_Tool.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "BRepTools.hxx"
#include "BRep_Builder.hxx"
#include "BRepTools_WireExplorer.hxx"
#include "BRepGProp.hxx"

#include "Interface_Static.hxx"

#include "MBTool.hpp"

#include "NCollection_IndexedDataMap.hxx"

#include "Poly_Array1OfTriangle.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_PolygonOnTriangulation.hxx"

#include "STEPControl_Reader.hxx"

#include "TColgp_Array1OfPnt.hxx"
#include "TColStd_Array1OfInteger.hxx"
#ifndef OCE
  #include "TopoDS_Shape.hxx"
#else
  #include "TopoDS_Shape.lxx"
#endif
#include "TopoDS.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Edge.hxx"
#include "TopLoc_Location.hxx"
#include "TopTools_HSequenceOfShape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"

typedef NCollection_IndexedDataMap<TopoDS_Face, moab::EntityHandle, TopTools_ShapeMapHasher> MapFaceToSurface;

typedef NCollection_IndexedDataMap<TopoDS_Edge, moab::EntityHandle, TopTools_ShapeMapHasher> MapEdgeToCurve;

struct TriangulationWithLocation {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation;
};

struct FacetingTolerance {
  float tolerance;
  bool is_relative;

  FacetingTolerance(float tol, bool is_absolute = false)
    : tolerance(tol), is_relative(!is_absolute) {}
};

void sew_and_facet(TopoDS_Shape &shape, const FacetingTolerance &facet_tol, MBTool &mbtool,
                   MaterialsMap &mat_map, std::string single_material = "");
void brep_faceter(std::string brep_file, std::string json_file,
                  const FacetingTolerance &facet_tol, std::string h5m_file, bool add_mat_ids);

#endif // BREP_FACETER_HH
