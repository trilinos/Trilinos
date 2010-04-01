
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
  Tifpack_Chebyshev.hpp
  Tifpack_Chebyshev_decl.hpp
  Tifpack_Chebyshev_def.hpp
  Tifpack_Condest.hpp
  Tifpack_CondestType.hpp
  Tifpack_ConfigDefs.hpp
  Tifpack_CreateOverlapGraph.hpp
  Tifpack_Diagonal.hpp
  Tifpack_Diagonal_decl.hpp
  Tifpack_Diagonal_def.hpp
  Tifpack_ExplicitInstantiationHelpers.hpp
  Tifpack_Factory.hpp
  Tifpack_Heap.hpp
  Tifpack_ILUT.hpp
  Tifpack_ILUT_decl.hpp
  Tifpack_ILUT_def.hpp
  Tifpack_IlukGraph.hpp
  Tifpack_OverlapGraph.hpp
  Tifpack_Parameters.hpp
  Tifpack_Preconditioner.hpp
  Tifpack_Relaxation.hpp
  Tifpack_Relaxation_decl.hpp
  Tifpack_Relaxation_def.hpp
  Tifpack_RILUK.hpp
  Tifpack_RILUK_decl.hpp
  Tifpack_RILUK_def.hpp
  Tifpack_ScalingType.hpp
  Tifpack_Version.hpp
  )

SET(SOURCES ${SOURCES}
  Tifpack_Chebyshev.cpp
  Tifpack_Diagonal.cpp
  Tifpack_Factory.cpp
  Tifpack_ILUT.cpp
  Tifpack_Parameters.cpp
  Tifpack_Relaxation.cpp
  Tifpack_RILUK.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  tifpack
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

