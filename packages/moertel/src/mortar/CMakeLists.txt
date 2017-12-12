

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

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

#
# Core files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  Moertel_Tolerances.hpp
  mrtr_pnode.H
  mrtr_function.H
  mrtr_point.H
  mrtr_functions.H
  mrtr_lm_selector.H
  mrtr_projector.H
  mrtr_segment_bilinearquad.H
  mrtr_segment_bilineartri.H
  mrtr_segment.H
  mrtr_segment_linear1D.H
  mrtr_node.H
  mrtr_overlap.hpp
  mrtr_overlap_Def.hpp
  mrtr_convexhull_Def.hpp
  mrtr_overlap_utils_Def.hpp
  mrtr_utils.H
  )

APPEND_SET(SOURCES
  mrtr_point.cpp
  mrtr_functions.cpp
  mrtr_projector3D.cpp
  mrtr_projector.cpp
  mrtr_interface_complete.cpp
  mrtr_segment_bilinearquad.cpp
  mrtr_interface_integrate3D.cpp
  mrtr_segment_bilineartri.cpp
  mrtr_segment.cpp
  mrtr_segment_linear1D.cpp
  mrtr_node.cpp
  mrtr_overlap.cpp
  mrtr_pnode.cpp
  mrtr_utils.cpp
  )

# Epetra specific

IF (HAVE_MOERTEL_EPETRA)

APPEND_SET(HEADERS
  mrtr_integrator.H
  mrtr_interface.H
  mrtr_manager.H
  mrtr_solver.H
  mrtr_ml_preconditioner.H
  )

APPEND_SET(SOURCES
  mrtr_integrator.cpp
  mrtr_interface_complete.cpp
  mrtr_interface_integrate3D.cpp
  mrtr_interface_integrate.cpp
  mrtr_interface_project.cpp
  mrtr_interface_tools.cpp
  mrtr_manager.cpp
  mrtr_manager_solve.cpp
  mrtr_solver.cpp
  mrtr_ml_preconditioner.cpp
  )

ENDIF ()

IF (HAVE_MOERTEL_TPETRA)

	APPEND_SET(HEADERS
		Moertel_ExplicitTemplateInstantiation.hpp
		Moertel_IntegratorT.hpp
		Moertel_IntegratorT_Def.hpp
		Moertel_InterfaceT.hpp
		Moertel_InterfaceT_Complete_Def.hpp
		Moertel_InterfaceT_Project_Def.hpp
		Moertel_InterfaceT_Integrate_Def.hpp
		Moertel_InterfaceT_Integrate3D_Def.hpp
		Moertel_InterfaceT_Tools_Def.hpp
		Moertel_ManagerT.hpp
		Moertel_ManagerT_Def.hpp
		Moertel_UtilsT.hpp
		Moertel_UtilsT_Def.hpp
	)
	
	APPEND_SET(SOURCES
		Moertel_ManagerT.cpp
		Moertel_InterfaceT.cpp
		Moertel_IntegratorT.cpp
	)

ENDIF ()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  moertel
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
