tribits_subpackage(StaticOnly)

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(
  static_only_lib
  STATIC
  HEADERS StaticOnly.hpp
  SOURCES StaticOnly.cpp
  )
  # NOTE: The library static_only_lib will be lined against the library
  # shared_only_lib from the upstream subpackage SharedOnly.  It is just fine
  # to link a static library against an upstream shared library but the
  # opposite is not true.

tribits_subpackage_postprocess()
