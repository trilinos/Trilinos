#
# A) Define the package
#

TRIBITS_PACKAGE(Gtest)

#
# B) Set up package-specific options
#

# 
# Disable Gtest's use of pthread.
#
ADD_DEFINITIONS(-DGTEST_HAS_PTHREAD=0)

#
# C) Add the libraries, tests, and examples
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS 
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest/gtest.h
  )

APPEND_SET(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/gtest/gtest-all.cc
  )

TRIBITS_ADD_LIBRARY(
  gtest
  NOINSTALLHEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#install gtest header into the correct subdir.
INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR}/gtest/)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
