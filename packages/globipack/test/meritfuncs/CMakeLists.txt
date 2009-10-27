
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  GlobiPack_TestLagrPolyMeritFunc1D.hpp
  GlobiPack_TestLagrPolyMeritFunc1D_decl.hpp
  GlobiPack_TestLagrPolyMeritFunc1D_def.hpp
  )

SET(SOURCES ${SOURCES}
  GlobiPack_TestLagrPolyMeritFunc1D.cpp
  )

PACKAGE_ADD_LIBRARY(
  globipack_test_meritfuncs
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  TESTONLY
  DEPLIBS globipack
  )

ADD_SUBDIRECTORY(UnitTests)
