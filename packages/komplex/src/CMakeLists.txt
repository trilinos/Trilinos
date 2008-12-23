
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

#SET(HAVE_BOOL ON)
#SET(HAVE_STL ON)

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

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# Core files
#

SET(HEADERS ${HEADERS}
  Komplex_LinearProblem.h
  azk_komplex.h
  Komplex_Version.h
  )

SET(SOURCES ${SOURCES}
Komplex_LinearProblem.cpp
  azk_create_linsys.c
  azk_create_precon.c
  azk_destroy_linsys.c
  azk_destroy_precon.c
  azk_create_matrix.c
  azk_create_vector.c
  azk_destroy_matrix.c
  azk_destroy_vector.c 
  )

#
# Comp files
#

SET(SOURCES ${SOURCES}
  azk_extract_solution.c
  azk_matvec.c
  azk_precon.c
  azk_permute_ri.c
  azk_residual_norm.c 
  )


#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  komplex
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
