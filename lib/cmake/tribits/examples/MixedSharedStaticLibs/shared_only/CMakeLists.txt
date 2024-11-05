tribits_subpackage(SharedOnly)

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(
  shared_only_lib
  SHARED
  HEADERS SharedOnly.hpp
  SOURCES SharedOnly.cpp
  )

tribits_add_test_directories(test)

tribits_subpackage_postprocess()
