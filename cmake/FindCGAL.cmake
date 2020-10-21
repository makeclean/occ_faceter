find_path(CGAL_CMAKE_CONFIG NAMES CGALConfig.cmake
  HINTS ${CGAL_ROOT}
  PATHS ENV LD_LIBRARY_PATH
  # ubuntu path
  PATHS /usr/lib/x86_64-linux-gnu/
  PATHS /usr/local/Cellar/cgal/4.12/lib/
  PATHS /usr/include/CGAL
  PATH_SUFFIXES lib Lib cmake cmake/CGAL
  NO_DEFAULT_PATH)

if(CGAL_CMAKE_CONFIG_STREQUAL "CGAL_CMAKE_CONFIG-NOTFOUND")
  set(CGAL_FOUND FALSE)
else()
  set(CGAL_FOUND TRUE)
  message(STATUS "Found CGAL in ${CGAL_CMAKE_CONFIG}")
  include(${CGAL_CMAKE_CONFIG}/CGALConfig.cmake)
endif()
