tribits_subpackage(C)

tribits_add_executable( c_util
  NOEXEPREFIX  NOEXESUFFIX
  SOURCES  c_util.cpp
  )

include(${CMAKE_CURRENT_LIST_DIR}/ShowLibErrors.cmake)

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
tribits_add_library(pws_c
  SOURCES C.cpp
  HEADERS wsp_c/C.hpp
  HEADERS_INSTALL_SUBDIR  wsp_c 
  ${TAL_EXTRALIB_ARGS} # Used to show warnings and errors
  )

tribits_add_test_directories(tests)

tribits_subpackage_postprocess()
