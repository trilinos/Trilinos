#
# A) Define the subpackage
#
tribits_subpackage(A)

#
# B) Set up subpackage-specific options
#
# Typically there are none or are few as most options are picked up from the
# parent package's CMakeLists.txt file!

#
# C) Add the libraries, tests, and examples
#

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
tribits_add_library(pws_a
  SOURCES A.cpp
  HEADERS A.hpp
  NOINSTALLHEADERS 
)

tribits_add_test_directories(tests)

#
# D) Do standard post processing
#
tribits_subpackage_postprocess()
