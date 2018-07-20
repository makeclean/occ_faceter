find_path(CGAL_CMAKE_CONFIG NAMES CGALConfig.cmake
  HINTS ${CGAL_ROOT}
  PATHS ENV LD_LIBRARY_PATH
  # ubuntu path
  PATHS /usr/lib/x86_64-linux-gnu/
  PATH_SUFFIXES lib Lib cmake cmake/CGAL
  NO_DEFAULT_PATH)

message(STATUS "Found CGAL in ${CGAL_CMAKE_CONFIG}")

include(${CGAL_CMAKE_CONFIG}/CGALConfig.cmake)
