TRIBITS_SUBPACKAGE(StaticOnly)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_LIBRARY(
  static_only_lib
  STATIC
  HEADERS StaticOnly.hpp
  SOURCES StaticOnly.cpp
  )
  # NOTE: The library static_only_lib will be lined against the libray
  # shared_only_lib from the upstream subpackage SharedOnly.  It is just fine
  # to link a static library against an upstream shared library but the
  # opposite is not true.

TRIBITS_SUBPACKAGE_POSTPROCESS()
