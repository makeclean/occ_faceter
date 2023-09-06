
IF(BUILD_ALL_TOOLS)
  ADD_EXECUTABLE(vertex_insert_test test_vertex_inserter.cc)
  TARGET_LINK_LIBRARIES(vertex_insert_test vertex_inserter ${MOAB_LIBRARIES})
  ADD_TEST(NAME vertex_insert_test COMMAND vertex_insert_test)
ENDIF(BUILD_ALL_TOOLS)

ADD_EXECUTABLE(faceter_test test_faceter.cc)
TARGET_LINK_LIBRARIES(faceter_test brep_faceter)
CONFIGURE_FILE(gluedCompSolid.brep . COPYONLY)
CONFIGURE_FILE(gluedCompSolid_materials.txt . COPYONLY)
CONFIGURE_FILE(dagmc_canon.py . COPYONLY)
CONFIGURE_FILE(test_occ_faceter_with_cubit_output.sh . COPYONLY)
CONFIGURE_FILE(cubit_output.h5m . COPYONLY)
ADD_TEST(NAME faceter_test COMMAND bash test_occ_faceter_with_cubit_output.sh)

ADD_EXECUTABLE(cylinder_faceting_test test_cylinder_faceting.cc)
TARGET_LINK_LIBRARIES(cylinder_faceting_test brep_faceter ${OCE_LINK_LIBS})
ADD_TEST(NAME cylinder_faceting_test COMMAND cylinder_faceting_test)
