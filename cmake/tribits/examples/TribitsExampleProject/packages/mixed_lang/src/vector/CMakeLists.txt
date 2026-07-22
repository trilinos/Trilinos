tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

tribits_add_library(
  mixedlang_vector
  HEADERONLY
  HEADERS Vector_Lite.hh
  HEADERS_INSTALL_SUBDIR vector
  )
