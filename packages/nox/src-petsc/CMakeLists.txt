
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base library
APPEND_SET(HEADERS
  NOX_Petsc.H
  NOX_Petsc_Vector.H
  NOX_Petsc_Group.H
  NOX_Petsc_SharedJacobian.H
  NOX_Petsc_Options.H
  NOX_Petsc_Interface.H
  )

APPEND_SET(SOURCES
  NOX_Petsc_Vector.C
  NOX_Petsc_Group.C
  NOX_Petsc_SharedJacobian.C
  NOX_Petsc_Options.C
  )

IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_PETSC)
  PACKAGE_ADD_LIBRARY(
    noxpetsc
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    DEPLIBS nox
    )
ENDIF()