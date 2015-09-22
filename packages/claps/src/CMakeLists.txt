
INCLUDE(CheckIncludeFileCXX)

#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS 
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core files
#


APPEND_SET(HEADERS
  Claps_ConfigDefs.hpp
  )

APPEND_SET(SOURCES
  Claps_Version.cpp
  )

IF (TPL_ENABLE_MPI)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/preconditioners)
  #
  # Extra files
  #

  APPEND_SET(HEADERS
    preconditioners/CLOP_constraint.hpp
    preconditioners/CLOP_graph.hpp
    preconditioners/CLOP_solver.hpp
    preconditioners/CLOP_sub.hpp
    preconditioners/CRD_utils.hpp
    preconditioners/EpetraExtCD_MatrixMatrix.hpp
    preconditioners/my_feti_sparse_solver.hpp
    preconditioners/myzero.hpp
    preconditioners/sparse_lu.hpp
    preconditioners/CLIP_solver.hpp
    preconditioners/sparse_lu2.hpp
    preconditioners/solver_crd.hpp
    preconditioners/CLIP_solver2.hpp
    preconditioners/CRS_serial.hpp
    preconditioners/preconditioner_crd.hpp
    preconditioners/krylov_crd.hpp
    )

  APPEND_SET(SOURCES
    preconditioners/CLOP_constraint.cpp
    preconditioners/CLOP_graph.cpp
    preconditioners/CLOP_solver.cpp
    preconditioners/CLOP_sub.cpp
    preconditioners/CRD_utils.cpp
    preconditioners/EpetraExtCD_MatrixMatrix.cpp
    preconditioners/sparse_lu.cpp
    preconditioners/CLIP_solver.cpp
    preconditioners/sparse_lu2.cpp
    preconditioners/CLIP_solver2.cpp
    preconditioners/CRS_serial.cpp
    preconditioners/krylov_crd.cpp
    )

  #
  # Extra Salinas files
  #
  IF (HAVE_SCLAPS)
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/salinas)
    
    APPEND_SET(HEADERS
      salinas/SCLOP_solver.hpp
      salinas/SCLIP_solver.hpp
      salinas/CLIP_interface.H
      salinas/CLOP_interface.H
      salinas/MpcLocal.H
      )

    APPEND_SET(SOURCES
      salinas/SCLOP_interface.cpp
      salinas/SCLOP_solver.cpp
      salinas/SCLIP_interface.cpp
      salinas/SCLIP_solver.cpp
      )

  ENDIF()


ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  claps
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
