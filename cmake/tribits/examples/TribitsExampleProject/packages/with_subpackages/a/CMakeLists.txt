#
# A) Define the subpackage
#
TRIBITS_SUBPACKAGE(A)

#
# B) Set up subpackage-specific options
#
# Typically there are none or are few as most options are picked up from the
# parent package's CMakeLists.txt file!

#
# C) Add the libraries, tests, and examples
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(pws_a
  SOURCES A.cpp
  HEADERS A.hpp
  NOINSTALLHEADERS 
)

TRIBITS_ADD_TEST_DIRECTORIES(tests)

#
# D) Do standard postprocessing
#
TRIBITS_SUBPACKAGE_POSTPROCESS()
