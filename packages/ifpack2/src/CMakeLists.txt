
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
  Ifpack2_Chebyshev.hpp
  Ifpack2_Chebyshev_decl.hpp
  Ifpack2_Chebyshev_def.hpp
  Ifpack2_Condest.hpp
  Ifpack2_CondestType.hpp
  Ifpack2_ConfigDefs.hpp
  Ifpack2_CreateOverlapGraph.hpp
  Ifpack2_Diagonal.hpp
  Ifpack2_Diagonal_decl.hpp
  Ifpack2_Diagonal_def.hpp
  Ifpack2_ExplicitInstantiationHelpers.hpp
  Ifpack2_Factory.hpp
  Ifpack2_Heap.hpp
  Ifpack2_ILUT.hpp
  Ifpack2_ILUT_decl.hpp
  Ifpack2_ILUT_def.hpp
  Ifpack2_IlukGraph.hpp
  Ifpack2_OverlapGraph.hpp
  Ifpack2_Parameters.hpp
  Ifpack2_Preconditioner.hpp
  Ifpack2_Relaxation.hpp
  Ifpack2_Relaxation_decl.hpp
  Ifpack2_Relaxation_def.hpp
  Ifpack2_RILUK.hpp
  Ifpack2_RILUK_decl.hpp
  Ifpack2_RILUK_def.hpp
  Ifpack2_ScalingType.hpp
  Ifpack2_BorderedOperator.hpp
  Ifpack2_Version.hpp
  )

SET(SOURCES ${SOURCES}
  Ifpack2_Chebyshev.cpp
  Ifpack2_Diagonal.cpp
  Ifpack2_Factory.cpp
  Ifpack2_ILUT.cpp
  Ifpack2_Parameters.cpp
  Ifpack2_Relaxation.cpp
  Ifpack2_RILUK.cpp
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

