
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

# TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

#
# All headers
#

#
# core headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/core)
APPEND_SET(HEADERS
  core/ROL_Algorithms.hpp
  )

#
# vspace headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vspace)
APPEND_SET(HEADERS
  vspace/ROL_stl.hpp
  vspace/ROL_teuchosarray.hpp
  )

#
# vector headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
APPEND_SET(HEADERS
  vector/ROL_Vector.hpp
  vector/ROL_StdVector.hpp
  vector/ROL_EpetraMultiVector.hpp
  )

#
# function headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
APPEND_SET(HEADERS
  function/ROL_Objective.hpp
  function/ROL_ObjectiveDef.hpp
)

#
# step headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/step)
APPEND_SET(HEADERS
  step/ROL_Step.hpp
  step/ROL_TrustRegionStep.hpp
  step/ROL_TrustRegion.hpp
  step/ROL_LineSearchStep.hpp
  step/ROL_LineSearch.hpp
  step/ROL_Secant.hpp
  step/ROL_lBFGS.hpp
  step/ROL_lDFP.hpp
  step/ROL_lSR1.hpp
  step/ROL_BarzilaiBorwein.hpp
  step/ROL_Krylov.hpp
  )

#
# algorithm headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithm)
APPEND_SET(HEADERS
  algorithm/ROL_Algorithm.hpp
  )

#
# status test headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/status)
APPEND_SET(HEADERS
  status/ROL_StatusTest.hpp
  )

#
# zoo headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/zoo)
APPEND_SET(HEADERS
  zoo/ROL_Types.hpp
  zoo/ROL_TestObjectives.hpp
  zoo/ROL_Rosenbrock.hpp
  zoo/ROL_FreudensteinRoth.hpp
  zoo/ROL_Beale.hpp
  zoo/ROL_Powell.hpp
  zoo/ROL_SumOfSquares.hpp
  zoo/ROL_LeastSquares.hpp
  zoo/ROL_PoissonControl.hpp
  )


#
# All sources
#

APPEND_SET(SOURCES
  ROL_Version.cpp
  )


IF (TPL_ENABLE_MPI)
ENDIF()


#
# C) Define the targets for package's library(s)
# 

  TRIBITS_ADD_LIBRARY(
    rol
    HEADERS ${HEADERS}
    SOURCES ${SOURCES}
    )
