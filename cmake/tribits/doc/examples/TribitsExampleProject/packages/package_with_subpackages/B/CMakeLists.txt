#
# A) Define the subpackage
#
TRIBITS_SUBPACKAGE(SubpackageB)

#
# B) Set up package-specific options
#
# None!

#
# C) Add the libraries, tests, and examples
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(pws_b
  SOURCES B.cpp
  HEADERS B.hpp
  NOINSTALLHEADERS 
)

TRIBITS_ADD_TEST_DIRECTORIES(tests)

#
# D) Do standard postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
