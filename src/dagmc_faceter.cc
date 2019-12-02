#include <iostream>
#include <array>

#include "BRep_Tool.hxx"
#include "BRepMesh_IncrementalMesh.hxx"
#include "Poly_Array1OfTriangle.hxx"
#include "TColgp_Array1OfPnt.hxx"
#include "TColStd_Array1OfInteger.hxx"
#include "Poly_Triangulation.hxx"
#include "Poly_PolygonOnTriangulation.hxx"

#include "Interface_Static.hxx"
//#include "IFSelect_PrintCount.hxx"
#include "STEPControl_Reader.hxx"

#include "TopoDS.hxx"
#include "TopoDS_Shape.hxx"
#include "TopoDS_Face.hxx"
#include "TopoDS_Edge.hxx"

#include "TopTools_HSequenceOfShape.hxx"
#include "TopTools_IndexedMapOfShape.hxx"
#include "TopExp_Explorer.hxx"
#include "TopExp.hxx"

//#include "TopLoc_TrsfPtr.hxx"
#include "TopLoc_Location.hxx"

#include "BRepOffsetAPI_Sewing.hxx"

#include "XSControl_WorkSession.hxx"
#include "Interface_InterfaceModel.hxx"
#include "StepRepr_Representation.hxx"
#include "TCollection_HAsciiString.hxx"

#include "MBTool.hpp"

MBTool *mbtool = new MBTool();

float facet_tol = 0.;

struct FaceterData {
  TopLoc_Location loc;
  Handle(Poly_Triangulation) triangulation;
};

//
facet_data get_facets(TopoDS_Face currentFace) {
  TopLoc_Location loc;
  
  //  facets.SetControlSurfaceDeflection(false);
   BRepMesh_IncrementalMesh facets(currentFace,facet_tol,false,0.5);
    //facets.theAngDeflection(0.5);
    //facets.theShape(currentFace);
    //facets.isRelative(false);
    //facets.theLinDeflection(1e-3f);
   facets.Perform();
    
   facet_data facetData;
  
  Handle(Poly_Triangulation) triangles = BRep_Tool::Triangulation(currentFace,loc);
  if(triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return facetData;
  } else {
    // retrieve facet data

    gp_Trsf local_transform = loc;
    const TColgp_Array1OfPnt &nodes = triangles->Nodes();
    for (int i = nodes.Lower() ; i <= nodes.Upper() ; i++) {
	Standard_Real x,y,z;
	nodes(i).Coord(x,y,z);
	local_transform.Transforms(x,y,z);
	std::array<double,3> coordinates;
	coordinates[0]=x,coordinates[1]=y,coordinates[2]=z;
	facetData.coords.push_back(coordinates);
	//std::cout << i << " " << x << " " << y << " " << z << std::endl;
    }
    // copy the facet_data
    std::array<int,3> conn;
    const Poly_Array1OfTriangle &tris = triangles->Triangles();
    std::cout << "Face has " << tris.Length() << " triangles" <<  std::endl;
    for (int i = tris.Lower() ; i <= tris.Upper() ; i++) {
	// get the node indexes for this triangle
	Poly_Triangle tri = tris(i);

	// reverse the triangle orientation if the face is reversed
	if (currentFace.Orientation() != TopAbs_FORWARD)
	  tri.Get(conn[2], conn[1], conn[0]);
	else
	  tri.Get(conn[0], conn[1], conn[2]);

	facetData.connectivity.push_back(conn);
    }	 
  return facetData;
  }
}

facet_data get_facets(TopoDS_Face currentFace, TopoDS_Edge currentEdge) {
  TopLoc_Location loc;
  //BRepMesh_IncrementalMesh facets;
  //  facets.SetControlSurfaceDeflection(false);
  #ifdef OCE_BUILD
    facets.SetAngle(0.5);
    facets.SetShape(currentFace);
    facets.SetRelative(false);
    facets.SetDeflection(1e-4f);
  #endif
  #ifdef OCC_BUILD
    BRepMesh_IncrementalMesh facets(currentFace,facet_tol,false,0.5);
    //facets.theAngDeflection(0.5);
    //facets.theShape(currentFace);
    //facets.isRelative(false);
    //facets.theLinDeflection(1e-3f);
  #endif
 
  facets.Perform();


  facet_data facetData;
  
  Handle(Poly_Triangulation) triangles = BRep_Tool::Triangulation(currentFace,loc);
  Handle(Poly_PolygonOnTriangulation) edges = BRep_Tool::PolygonOnTriangulation(currentEdge,triangles,loc);
  return facetData;
}


// get all the curves of a given shape
void get_edges(TopoDS_Shape shape) {
  TopTools_IndexedMapOfShape edges;
  TopExp::MapShapes(shape,TopAbs_EDGE,edges);
  for ( int i = 1 ; i<=edges.Extent() ; i++ ) {
    TopoDS_Edge currentEdge = TopoDS::Edge(edges(i));
    //    facet_data facets = get_facets(currentEdge);
    //   std::cout << i << std::endl;    
  }
  //for ( TopExp_Explorer ex(shape, TopAbs_EDGE) ; ex.More() ; ex.Next() ) {
  //  TopoDS_Edge currentEdge = TopoDS::Edge(ex.Current());
  //    facet_data facets = get_facets(currentEdge);
  return;
}

// get the triangulation for the current face
FaceterData get_triangulation(TopoDS_Face currentFace) {
  TopLoc_Location loc;
  //BRepMesh_IncrementalMesh facets;
  //  facets.SetControlSurfaceDeflection(false);
  #ifdef OCE_BUILD
    facets.SetAngle(0.5);
    facets.SetShape(currentFace);
  facets.SetRelative(false);
  facets.SetDeflection(facet_tol);
  #endif
  #ifdef OCC_BUILD
    BRepMesh_IncrementalMesh facets(currentFace,facet_tol,false,0.5);
    //facets.theAngDeflection(0.5);
    //facets.theShape(currentFace);
    //facets.isRelative(false);
    //facets.theLinDeflection(1e-3f);
  #endif
  facets.Perform();

  Handle(Poly_Triangulation) triangles = BRep_Tool::Triangulation(currentFace,loc);
  FaceterData data;
  data.loc = loc;
  data.triangulation = triangles;
  return data;
}

facet_data make_surface_facets(TopoDS_Face currentFace, FaceterData facetData) {
  facet_data facets_for_moab;

  Handle(Poly_Triangulation) triangles = facetData.triangulation;
  TopLoc_Location loc = facetData.loc;  
  
  if(triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return facets_for_moab;
  } else {
    // retrieve facet data

    gp_Trsf local_transform = loc;
    const TColgp_Array1OfPnt &nodes = triangles->Nodes();
    for (int i = nodes.Lower() ; i <= nodes.Upper() ; i++) {
	Standard_Real x,y,z;
	nodes(i).Coord(x,y,z);
	local_transform.Transforms(x,y,z);
	std::array<double,3> coordinates;
	coordinates[0]=x,coordinates[1]=y,coordinates[2]=z;
	facets_for_moab.coords.push_back(coordinates);
	//std::cout << i << " " << x << " " << y << " " << z << std::endl;
    }
    // copy the facet_data
    std::array<int,3> conn;
    const Poly_Array1OfTriangle &tris = triangles->Triangles();
    std::cout << "Face has " << tris.Length() << " triangles" <<  std::endl;
    for (int i = tris.Lower() ; i <= tris.Upper() ; i++) {
	// get the node indexes for this triangle
	Poly_Triangle tri = tris(i);

	// reverse the triangle orientation if the face is reversed
	if (currentFace.Orientation() != TopAbs_FORWARD)
	  tri.Get(conn[2], conn[1], conn[0]);
	else
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
  Handle(Poly_PolygonOnTriangulation) edges =  BRep_Tool::PolygonOnTriangulation(currentEdge,triangles,loc);
  
  if(triangles.IsNull()) {
    std::cout << "No facets for surface" << std::endl;
    return edges_for_moab;
  } else {
    // retrieve facet data
    std::vector<int> conn;
    //TColStd_Array1OfInteger
    const TColStd_Array1OfInteger &lines = edges->Nodes();
    //std::cout << "Edge has " << lines.Length() << " length" <<  std::endl;
    for (int i = lines.Lower() ; i <= lines.Upper() ; i++) {
      conn.push_back(lines(i));
    }
    edges_for_moab.connectivity = conn;
    return edges_for_moab;
  }
}

// for a given shape - get the faces and edges
// then facet the face and also extract edges
void get_facets_for_shape(TopoDS_Shape shape) {
  int j = 0;
  moab::EntityHandle vol;

  //#pragma omp critical
  moab::ErrorCode rval = mbtool->make_new_volume(vol);
  
  for( TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next() ) {
    TopoDS_Face currentFace = TopoDS::Face( ex.Current());
    //facet_data facets = get_facets(currentFace);
    // get the triangulation for the current face
    FaceterData data = get_triangulation(currentFace);
    // make facets for current face
    facet_data facets = make_surface_facets(currentFace,data);
    
    TopTools_IndexedMapOfShape edges;
    TopExp::MapShapes(currentFace,TopAbs_EDGE,edges);
    std::vector<edge_data> edge_collection;
    for ( int i = 1 ; i<=edges.Extent() ; i++ ) {
      TopoDS_Edge currentEdge = TopoDS::Edge(edges(i));
      // make the edge facets
      edge_data edges = make_edge_facets(currentFace,currentEdge,data);
      edge_collection.push_back(edges);
      //    std::cout << j << " " <<  i << std::endl;
    }
    //std::cout << edge_collection.size() << std::endl;
    //mbtool->add_surface(vol,facets);//,edge_collection);
    //#pragma omp critical
    mbtool->add_surface(vol,facets,edge_collection);
    //std::cout << j << " " << edges.Extent() << std::endl;
  }
  return;
}

// facet all the volumes
void facet_all_volumes(Handle_TopTools_HSequenceOfShape shape_list){
  int count = shape_list->Length();
  //#pragma omp parallel for 
  for ( int i = 1 ; i <= count ; i++ ) {
    TopoDS_Shape shape = shape_list->Value(i);
    get_facets_for_shape(shape); // get the edges and 
    //get_faces(shape);
    //get_edges(shape);
  }
  return;
}

STEPControl_Reader *step;

//dump all labels
void dump_labels() {
  const Handle(XSControl_WorkSession)& theSession = step->WS();
  const Handle(Interface_InterfaceModel)& theModel = theSession->Model();

  Standard_Integer nb = theModel->NbEntities();
  for(Standard_Integer i=1; i<=nb; i++)
    {

      Handle(StepRepr_Representation) ent =
	Handle(StepRepr_Representation)::DownCast
	(theModel->Value(i));
      if (ent.IsNull()) continue;
      if (ent->Name().IsNull()) continue;
      std::cout << ent->Name()->ToCString() << std::endl;
    }
  return;
}

void sew_and_append(TopoDS_Shape shape, Handle(TopTools_HSequenceOfShape) &shape_list) {
  std::cout << "shape type: " << shape.ShapeType() << std::endl;
  // sew together all the curves
  BRepOffsetAPI_Sewing(1.0e-06, Standard_True);
  BRepOffsetAPI_Sewing sew;
  sew.Add(shape);
  sew.Perform();
  shape = sew.SewedShape();
  shape_list->Append(shape);
  return;
}


int main (int argc, char* argv[]) {

  std::string cad_file(argv[1]);
  std::string ftol(argv[2]);
  facet_tol = std::stod(ftol);
  std::string filename(argv[3]);

  moab::ErrorCode rval = mbtool->set_tags();
  
  //Interface_Static::SetIVal("read.step.product.mode",0);
  //Interface_Static::SetIVal("read.step.product.context",2);
  //Interface_Static::SetIVal("read.step.assembly.level",4);
  
  step = new STEPControl_Reader();  
  step->ReadFile(cad_file.c_str());

  //step->PrintCheckLoad(false,IFSelect_GeneralInfo);
  //step->PrintCheckLoad(false,IFSelect_ItemsByEntity);
  //step->PrintCheckLoad(false,IFSelect_CountByItem);
  step->PrintCheckLoad(false,IFSelect_ListByItem);
  step->ClearShapes();
  int count = step->NbRootsForTransfer();
  TopoDS_Shape shape;

  int r_count = step->TransferRoots();
  std::cout << "r_count: " << r_count << std::endl;
  std::cout << "n_shapes: " << step->NbShapes() << std::endl;
    
  Handle(TopTools_HSequenceOfShape) shape_list = new TopTools_HSequenceOfShape;
  std::cout << count << std::endl;
  for ( int i = 1 ; i <= count ; i++ ) {
    bool ok = step->TransferRoot(i);
    step->PrintCheckTransfer(false,IFSelect_CountByItem);
    if (ok) {
      shape = step->Shape(i);
      // if its a compound decend and get its children
      if ( shape.ShapeType() == 0 ) {
	TopoDS_Iterator it = TopoDS_Iterator(shape);
	while ( it.More() ) {
	  // if we are a volume
	  if ( it.Value().ShapeType() == 2 ) {
	    sew_and_append(it.Value(),shape_list);
	  } else {
	    std::cout << "Unknown shape type " << it.Value().ShapeType() << std::endl;
	  }
	  it.Next();
	}
      } else if ( shape.ShapeType() == 2) {
	// we are a normal volume insert into
	// the list
	sew_and_append(shape,shape_list);
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

