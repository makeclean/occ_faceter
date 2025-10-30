
#ADD_EXECUTABLE(vertex_insert_test test_vertex_inserter.cc)
#TARGET_LINK_LIBRARIES(vertex_insert_test vertex_inserter ${MOAB_LIBRARIES})
ADD_TEST(NAME vertex_insert_test COMMAND vertex_insert_test)

IF(OCE_FOUND OR OCC_FOUND )
ADD_EXECUTABLE(faceter_test test_faceter.cc)
TARGET_LINK_LIBRARIES(faceter_test brep_faceter)
CONFIGURE_FILE(gluedCompSolid.brep . COPYONLY)
CONFIGURE_FILE(gluedCompSolid_metadata.json . COPYONLY)
ADD_TEST(NAME faceter_test COMMAND faceter_test)
ENDIF(OCE_FOUND OR OCC_FOUND )

# mesh utils test
ADD_EXECUTABLE(surf_utils_test test_surf_utils.cc)
TARGET_LINK_LIBRARIES(surf_utils_test surface_mesh_utils)
CONFIGURE_FILE(test.h5m . COPYONLY)
CONFIGURE_FILE(test_2.h5m . COPYONLY)
ADD_TEST(NAME surf_utils_test COMMAND surf_utils_test)

