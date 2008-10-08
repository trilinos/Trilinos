
INCLUDE(TrilinosPackageLibraryHelpers)

#
# A) Package-specific configuration options
#

TRILINOS_PACKAGE_CONFIGURE_FILE(${PROJECT_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})


SET(CORE_H
  Trilinos_Util.h
  Trilinos_Util_ReadMatrixMarket2Epetra.h
  Trilinos_Util_CountMatrixMarket.h
  Trilinos_Util_CountTriples.h
  Trilinos_Util_CrsMatrixGallery.h
  Trilinos_Util_CommandLineParser.h
  )

SET(CORE
  Trilinos_Util_duscr_vbr.cpp
  Trilinos_Util_dusds_vbr.cpp
  Trilinos_Util_distrib_vbr_matrix.cpp
  Trilinos_Util_read_hb.cpp
  Trilinos_Util_svbrres.cpp
  Trilinos_Util_dusmm.cpp
  Trilinos_Util_scscmv.cpp
  Trilinos_Util_read_vec.cpp
  Trilinos_Util_write_vec.cpp
  Trilinos_Util_create_vbr.cpp
  Trilinos_Util_msr2vbr.cpp
  Trilinos_Util_scscres.cpp
  Trilinos_Util_distrib_msr_matrix.cpp
  Trilinos_Util_read_coo.cpp
  Trilinos_Util_smsrres.cpp
  Trilinos_Util_csrcsc.cpp
  Trilinos_Util_csrmsr.cpp
  Trilinos_Util_ssrcsr.cpp
  Trilinos_Util_coocsr.cpp
  Trilinos_Util_ReadMatrixMarket2Epetra.cpp
  Trilinos_Util_ReadTriples2Epetra.cpp
  Trilinos_Util_CountMatrixMarket.cpp
  Trilinos_Util_CountTriples.cpp
  Trilinos_Util_ReadHb2Epetra.cpp
  Trilinos_Util_ReadHpc2Epetra.cpp
  Trilinos_Util_ReadHb2EpetraVbr.cpp
  Trilinos_Util_GenerateCrsProblem.cpp
  Trilinos_Util_GenerateVbrProblem.cpp
  Trilinos_Util_CommandLineParser.cpp
  Trilinos_Util_CrsMatrixGallery.cpp
  )

SET(UTIL_H
  iohb.h
  Trilinos_Util_Version.h
  )

SET(UTIL
  iohb.cpp
  )

#
# C) Define the targets for package's library(s)
#

TRILINOS_PACKAGE_ADD_LIBRARY(
  triutils
  HEADERS
    ${CORE_H}
    ${UTIL_H}
  SOURCES
    ${CORE}
    ${UTIL}
  )

#
# D) Export the dependancy variables of this package for other packages to use
#

TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES()
