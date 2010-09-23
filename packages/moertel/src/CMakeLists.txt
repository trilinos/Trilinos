
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
# Core files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  mrtr_pnode.H
  mrtr_function.H
  mrtr_point.H
  mrtr_functions.H
  mrtr_integrator.H
  mrtr_lm_selector.H
  mrtr_projector.H
  mrtr_interface.H
  mrtr_segment_bilinearquad.H
  mrtr_segment_bilineartri.H
  mrtr_segment.H
  mrtr_manager.H
  mrtr_segment_linear1D.H
  mrtr_solver.H
  mrtr_ml_preconditioner.H
  mrtr_node.H
  mrtr_overlap.H
  )

APPEND_SET(SOURCES
  mrtr_convexhull.cpp
  mrtr_point.cpp
  mrtr_functions.cpp
  mrtr_projector3D.cpp
  mrtr_integrator.cpp
  mrtr_projector.cpp
  mrtr_interface_complete.cpp
  mrtr_segment_bilinearquad.cpp
  mrtr_interface_integrate3D.cpp
  mrtr_segment_bilineartri.cpp
  mrtr_interface_integrate.cpp
  mrtr_interface_project.cpp
  mrtr_segment.cpp
  mrtr_interface_tools.cpp
  mrtr_manager.cpp
  mrtr_segment_linear1D.cpp
  mrtr_manager_solve.cpp
  mrtr_solver.cpp
  mrtr_ml_preconditioner.cpp
  mrtr_node.cpp
  mrtr_overlap.cpp
  mrtr_overlap_utils.cpp
  mrtr_pnode.cpp
  )

#
# Util files
#

APPEND_SET(HEADERS
  mrtr_utils.H
  )

APPEND_SET(SOURCES
  mrtr_utils.cpp
  )

IF (TPL_ENABLE_MPI)
ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  moertel
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
