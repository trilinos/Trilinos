
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)
INCLUDE(TribitsLibraryMacros)

#
# A) Define the package
#

TRIBITS_PACKAGE(Gtest)

#
# B) Set up package-specific options
#

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
  NO_INSTALL_HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

#install gtest header into the correct subdir.
INSTALL(FILES ${HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${Trilinos_INSTALL_INCLUDE_DIR}/gtest/)



#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
