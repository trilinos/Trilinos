#
# A) Define the subpackage
#
tribits_subpackage(A)

#
# B) Set up subpackage-specific options
#

set(${PACKAGE_NAME}_SPECIAL_VALUE  3  CACHE STRING "Integer special value")
tribits_pkg_export_cache_var(${PACKAGE_NAME}_SPECIAL_VALUE)

#
# C) Add the libraries, tests, and examples
#

tribits_configure_file(${PACKAGE_NAME}_config.h)

tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})
tribits_add_library(pws_a
  SOURCES A.cpp
  HEADERS A.hpp ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

tribits_add_test_directories(tests)

#
# D) Do standard post processing
#
tribits_subpackage_postprocess()
