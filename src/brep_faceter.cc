#include "brep_faceter.hh"

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
  std::vector<int> senses;
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
void perform_faceting(const TopoDS_Face &face, const FacetingTolerance& facet_tol) {
  // This constructor calls Perform()
  BRepMesh_IncrementalMesh facets(face, facet_tol.tolerance, facet_tol.is_relative, 0.5);
  // todo need to get the wires that belong to the face, and query their sense relative
  // to the face
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
  MapFaceToSurface surfaceMap; // map of OCCFace to Surface entities
  MapEdgeToCurve curveMap; // map of OCCEdge to Curve entities

  // list unique faces, create empty surfaces, and build surfaceMap
  // for the purpose of performing the faceting in parallel
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
      mbtool.make_new_surface_tags(surface);
      surfaceMap.Add(face, surface);

      // investigate all the curves that belong to the shape
      TopExp_Explorer Ex;
      TopTools_IndexedMapOfShape edgeMap;
      TopExp::MapShapes(face,TopAbs_EDGE,edgeMap);
      std::cout << surface << " " << edgeMap.Extent() << std::endl;
      for (int i = 1 ; i <= edgeMap.Extent() ; i++ ) {
	TopoDS_Edge aEdge = TopoDS::Edge(edgeMap(i));
	moab::EntityHandle curve;
	mbtool.make_new_curve_tags(curve);
	curveMap.Add(aEdge,curve);
	// get the vertices for the edge
	TopExp_Explorer verts(aEdge,TopAbs_VERTEX);
	TopTools_IndexedMapOfShape vertexMap;
	TopExp::MapShapes(aEdge,TopAbs_VERTEX,vertexMap);
	std::cout << curve << " " << vertexMap.Extent() << std::endl;
	if(!verts.More()) {
	  std::cout << "No vertices" << std::endl;
	} else {
	  verts.ReInit();
	}
  
	for ( TopExp_Explorer verts ; verts.More();verts.Next()) {
	  TopoDS_Vertex vert = TopoDS::Vertex(verts.Current());
	  gp_Pnt p = BRep_Tool::Pnt(vert);
	  std::array<double,3> coord = {double(p.X()),double(p.Y()),double(p.Z())};
	  moab::EntityHandle vertex,vertex_set;
	  moab::ErrorCode rval = mbtool.add_vertex(coord,vertex,vertex_set);
	  mbtool.add_vertex_to_curve(curve,vertex_set);
	}
      }	
      // curve topology now exists
    }
  }
  // surface entity sets now exist
  
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
    //    get_curve_data(face,data); // get the curve data
    mbtool.add_facets_and_curves_to_surface(surface, data.facets, data.edge_collection);
  }

  // build a list of volumes for each material
  std::map<std::string, std::vector<moab::EntityHandle>> material_volumes;

  // create volumes and add surfaces to parents
  for (int i = 1; i <= count; i++) {
    const TopoDS_Shape &shape = shape_list.Value(i);

    moab::EntityHandle vol;
    mbtool.make_new_volume_tags(vol);
    
    // loop over all the surfaces adding the relationship between surfaces and volumes
    for (TopExp_Explorer ex(shape, TopAbs_FACE); ex.More(); ex.Next()) {
      const TopoDS_Face &face = TopoDS::Face(ex.Current());
      moab::EntityHandle surface = surfaceMap.FindFromKey(face);
      int sense = face.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
      mbtool.add_surface_to_volume(surface, vol, sense);
      // now should set the sense of each curve wrt to the surface
      for(TopExp_Explorer wires(face,TopAbs_WIRE); wires.More(); wires.Next()) {
	// explore the wires
	const TopoDS_Wire &wire = TopoDS::Wire(wires.Current());
	//BRepTools_WireExplorer WireEx();
	//WireEx.Init(wire, face);
	//BRepTools_WireExplorer WireEx(wire);
	//for (WireEx.Begin() ; WireEx.More(); WireEx.Next()) {
	for (BRepTools_WireExplorer WireEx(wire,face); WireEx.More(); WireEx.Next()) {
	  const TopoDS_Edge &edge = TopoDS::Edge(WireEx.Current());
	  int curve_sense = edge.Orientation() == TopAbs_REVERSED ? moab::SENSE_REVERSE : moab::SENSE_FORWARD;
	  std ::cout << &edge << " " << sense << std::endl;
	  mbtool.add_curve_to_surface(curveMap.FindFromKey(edge), surface, curve_sense);
	  // Process current edge
	}
      }
    }

    // now loopover the surfaces of the shape - and query the curve
    // ownership and sense 
    
    // curve senses now set
    
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
  mbtool.set_tags();
  sew_and_facet(shape, facet_tol, mbtool, materials_map);

  if (add_mat_ids)
    mbtool.add_mat_ids();

  mbtool.write_geometry(h5m_file.c_str());
}
