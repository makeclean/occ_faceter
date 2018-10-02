# cgal includes
INCLUDE_DIRECTORIES(${MOAB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${OPENCASCADE_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(catch)

ADD_LIBRARY(vertex_inserter vertex_inserter.cc)

ADD_LIBRARY(mbtool MBTool.cpp)

# occ faceter
LINK_DIRECTORIES("/usr/lib/x86_64-linux-gnu/")
IF(OCE_FOUND)
  ADD_EXECUTABLE(occ_faceter dagmc_faceter.cc)
  SET(OCE_LINK_LIBS TKSTEP TKXSBase TKernel TKMath TKMesh TKTopAlgo TKBRep TKSTEPBase)
  TARGET_LINK_LIBRARIES(occ_faceter mbtool vertex_inserter ${MOAB_LIBRARIES} ${OCE_LINK_LIBS})
  INSTALL(TARGETS occ_faceter DESTINATION bin)
ENDIF(OCE_FOUND)

ADD_LIBRARY(dagmc_topology dagmc_topology.cc)
# dagmc merger
ADD_EXECUTABLE(dagmc_merge dagmc_merge.cc)
TARGET_LINK_LIBRARIES(dagmc_merge dagmc_topology ${MOAB_LIBRARIES})

# slcing library
ADD_LIBRARY(dagmc_slicer dagmc_slicer.cc)
TARGET_LINK_LIBRARIES(dagmc_slicer ${MOAB_LIBRARIES} ${CGAL_LIBRARIES} boost_system gmp)

# dagmc slicer
ADD_EXECUTABLE(dagmc_slice_tool dagmc_slice_tool.cc)
TARGET_LINK_LIBRARIES(dagmc_slice_tool dagmc_slicer ${MOAB_LIBRARIES} ${CGAL_LIBRARY} boost_system boost_program_options gmp)

# dagmc 2d mesher
ADD_EXECUTABLE(dagmc_2d_mesh dagmc_2d_mesh.cc)
TARGET_LINK_LIBRARIES(dagmc_2d_mesh dagmc_slicer vertex_inserter ${MOAB_LIBRARIES} ${CGAL_LIBRARY} boost_system boost_program_options gmp)

# dagmc_overlap
ADD_EXECUTABLE(dagmc_overlap preper.cc)
TARGET_LINK_LIBRARIES(dagmc_overlap vertex_inserter ${MOAB_LIBRARIES} ${CGAL_LIBRARY} boost_system boost_program_options gmp mpfr)

INSTALL(TARGETS dagmc_slicer DESTINATION bin)
INSTALL(TARGETS dagmc_merge DESTINATION bin)
INSTALL(TARGETS dagmc_overlap DESTINATION bin)
INSTALL(TARGETS dagmc_2d_mesh DESTINATION bin)

ADD_SUBDIRECTORY(test)
