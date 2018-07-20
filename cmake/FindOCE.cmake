find_path(OCE_CMAKE_CONFIG NAMES OCEConfig.cmake
  HINTS ${OCE_ROOT}
  PATHS ENV LD_LIBRARY_PATH
  # ubuntu path
  PATHS /usr/lib/x86_64-linux-gnu/oce-0.17
  PATH_SUFFIXES lib Lib cmake cmake/OCE
  NO_DEFAULT_PATH)

message(STATUS "Found OCE in ${OCE_CMAKE_CONFIG}")

include(${OCE_CMAKE_CONFIG}/OCEConfig.cmake)
