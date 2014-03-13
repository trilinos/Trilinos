#
# A) Define the subpackage
#

TRIBITS_SUBPACKAGE(SubpackageC ENABLE_SHADOWING_WARNINGS  CLEANED )

#
# B) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# C) Add the libraries, tests, and examples
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(pws_c
  SOURCES C.cpp
  HEADERS C.hpp
)

TRIBITS_ADD_TEST_DIRECTORIES(tests)

#
# D) Do standard postprocessing
#

TRIBITS_SUBPACKAGE_POSTPROCESS()
