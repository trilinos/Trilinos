
INCLUDE(PackageLibraryMacros)

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base library
APPEND_SET(HEADERS
  NOX_Belos_MultiVector.H
  NOX_Belos_JacobianOperator.H
  NOX_Belos_PreconditionOperator.H
  NOX_Belos_Group.H
  )

APPEND_SET(SOURCES
  NOX_Belos_MultiVector.C
  NOX_Belos_JacobianOperator.C
  NOX_Belos_PreconditionOperator.C
  NOX_Belos_Group.C
  )

IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_BELOS)
  PACKAGE_ADD_LIBRARY(
    noxbelos
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
ENDIF()
