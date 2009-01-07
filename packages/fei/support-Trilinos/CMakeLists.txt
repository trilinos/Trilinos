
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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../)

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  fei_Aztec_BlockMap.hpp
  fei_AztecDMSR_Matrix.hpp
  fei_AztecDVBR_Matrix.hpp
  fei_Aztec_LinSysCore.hpp
  fei_Aztec_Map.hpp
  fei_Aztec_Vector.hpp
  cfei_aztec.h
  fei_Factory_Trilinos.hpp
  fei_Include_Trilinos.hpp
  fei_LinProbMgr_EpetraBasic.hpp
  fei_MatrixTraits_Epetra.hpp
  fei_Solver_Amesos.hpp
  fei_Solver_AztecOO.hpp
  fei_Trilinos_Helpers.hpp
  fei_trilinos_macros.hpp
  fei_VectorTraits_Epetra.hpp
  )

APPEND_SET(SOURCES
  fei_Aztec_BlockMap.cpp
  fei_AztecDMSR_Matrix.cpp
  fei_AztecDVBR_Matrix.cpp
  fei_Aztec_LinSysCore.cpp
  fei_Aztec_Map.cpp
  fei_Aztec_Vector.cpp
  cfei_aztec.cpp
  fei_Factory_Trilinos.cpp
  fei_LinProbMgr_EpetraBasic.cpp
  fei_Solver_Amesos.cpp
  fei_Solver_AztecOO.cpp
  fei_Trilinos_Helpers.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fei_trilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
