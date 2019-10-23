
ADD_EXECUTABLE(vertex_insert_test test_vertex_inserter.cc)
TARGET_LINK_LIBRARIES(vertex_insert_test vertex_inserter ${MOAB_LIBRARIES})

ADD_TEST(NAME vertex_insert_test COMMAND vertex_insert_test)
