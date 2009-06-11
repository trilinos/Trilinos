
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

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Kokkos_ConfigDefs.hpp
  Kokkos_BaseSparseMultiply.hpp
  Kokkos_BaseSparseSolve.hpp
  Kokkos_CisMatrix.hpp
  Kokkos_DenseMultiVector.hpp
  Kokkos_DenseVector.hpp
  Kokkos_HbMatrix.hpp
  Kokkos_MultiVector.hpp
  Kokkos_PackedSparseMultiply.hpp
  Kokkos_Permutation.hpp
  Kokkos_SparseOperation.hpp
  Kokkos_Vector.hpp
  )

#
# Util files
#

APPEND_SET(HEADERS
  Kokkos_Flops.hpp
  Kokkos_Time.hpp
  Kokkos_Fortran_wrappers.hpp
  Kokkos_f77func.hpp
  Kokkos_Version.hpp
  )

APPEND_SET(SOURCES
  Kokkos_dcrsmv.f
  Kokkos_dcrsmm.f
  )

APPEND_SET(HEADERS
  Kokkos_CompObject.hpp
  )

######################################

APPEND_SET(HEADERS
  )

APPEND_SET(SOURCES
  )

######################################
IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  kokkos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
