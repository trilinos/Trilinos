
#INCLUDE(TribitsLibraryMacros)

#SET(HEADERS "")
#SET(SOURCES "")

#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Base library
#APPEND_SET(HEADERS
#  NOX_Thyra.H
#  NOX_Thyra_Vector.H
#  NOX_Thyra_MultiVector.H
#  NOX_Thyra_Group.H
#  Thyra_NonlinearSolver_NOX.hpp
#  NOX_PrePostOperator_RowSumScaling.H
#  NOX_MeritFunction_Weighted.hpp
#  NOX_Solver_PseudoTransient.hpp
#  )

#APPEND_SET(SOURCES
#  NOX_Thyra_Vector.C
#  NOX_Thyra_MultiVector.C
#  NOX_Thyra_Group.C
#  Thyra_NonlinearSolver_NOX.cpp
#  NOX_PrePostOperator_RowSumScaling.C
#  NOX_MeritFunction_Weighted.cpp
#  NOX_Solver_PseudoTransient.cpp
#  )

#IF (NOX_ENABLE_ABSTRACT_IMPLEMENTATION_THYRA)
#  TRIBITS_ADD_LIBRARY(
#    noxthyra
#    HEADERS ${HEADERS}
#    SOURCES ${SOURCES}
#    DEPLIBS nox
#    )
#ENDIF()