
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base library
APPEND_SET(HEADERS
  NOX_Thyra.H
  NOX_Thyra_Vector.H
  NOX_Thyra_MultiVector.H
  NOX_Thyra_Group.H
  Thyra_NonlinearSolver_NOX.hpp
  )

APPEND_SET(SOURCES
  NOX_Thyra_Vector.C
  NOX_Thyra_MultiVector.C
  NOX_Thyra_Group.C
  Thyra_NonlinearSolver_NOX.cpp
  )

IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA)
  PACKAGE_ADD_LIBRARY(
    noxthyra
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()