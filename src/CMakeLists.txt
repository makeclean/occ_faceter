# cgal includes
INCLUDE_DIRECTORIES(${MOAB_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${OCE_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(catch)

ADD_LIBRARY(vertex_inserter vertex_inserter.cc)

ADD_LIBRARY(mbtool MBTool.cpp)

# occ faceter
LINK_DIRECTORIES("/usr/lib/x86_64-linux-gnu/")
ADD_EXECUTABLE(occ_faceter dagmc_faceter.cc)
SET(OCE_LINK_LIBS TKSTEP TKXSBase TKernel TKMath TKMesh TKTopAlgo TKBRep TKSTEPBase)
TARGET_LINK_LIBRARIES(occ_faceter mbtool vertex_inserter ${MOAB_LIBRARIES} ${OCE_LINK_LIBS})

ADD_LIBRARY(dagmc_topology dagmc_topology.cc)
# dagmc merger
ADD_EXECUTABLE(dagmc_merge dagmc_merge.cc)
TARGET_LINK_LIBRARIES(dagmc_merge dagmc_topology ${MOAB_LIBRARIES})

# dagmc slicer
ADD_EXECUTABLE(dagmc_slicer dagmc_slicer.cc)
TARGET_LINK_LIBRARIES(dagmc_slicer ${MOAB_LIBRARIES} ${CGAL_LIBRARY} boost_system boost_program_options gmp)

INSTALL(TARGETS dagmc_slicer DESTINATION bin)
INSTALL(TARGETS dagmc_merge DESTINATION bin)
INSTALL(TARGETS occ_faceter DESTINATION bin)

ADD_SUBDIRECTORY(test)
