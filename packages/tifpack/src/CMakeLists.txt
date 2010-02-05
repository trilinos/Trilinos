
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS)
SET(SOURCES)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Tifpack_ConfigDefs.hpp
  Tifpack_Chebyshev.hpp
  Tifpack_Condest.hpp
  Tifpack_CondestType.hpp
  Tifpack_ExplicitInstantiationHelpers.hpp
  Tifpack_Heap.hpp
  Tifpack_ILUT.hpp
  Tifpack_ILUT_decl.hpp
  Tifpack_ILUT_def.hpp
  Tifpack_IlukGraph.hpp
  Tifpack_Parameters.hpp
  Tifpack_Preconditioner.hpp
  Tifpack_Relaxation.hpp
  Tifpack_ScalingType.hpp
  Tifpack_Version.hpp
  )

SET(SOURCES ${SOURCES}
  Tifpack_Factory.cpp
  Tifpack_ILUT.cpp
  Tifpack_Parameters.cpp
  Tifpack_Relaxation.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tifpack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

