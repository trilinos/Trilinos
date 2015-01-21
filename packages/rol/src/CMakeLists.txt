
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
# vector headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/vector)
APPEND_SET(HEADERS
  vector/ROL_Vector.hpp
  vector/ROL_Vector_SimOpt.hpp
  vector/ROL_StdVector.hpp
  )

#
# function headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/function)
APPEND_SET(HEADERS
  function/ROL_Objective.hpp
  function/ROL_ObjectiveDef.hpp
  function/ROL_Objective_SimOpt.hpp
  function/ROL_Reduced_Objective_SimOpt.hpp
  function/ROL_Constraints.hpp
  function/ROL_EqualityConstraint.hpp
  function/ROL_EqualityConstraintDef.hpp
  function/ROL_EqualityConstraint_SimOpt.hpp
  function/ROL_BoundConstraint.hpp
  function/ROL_StdBoundConstraint.hpp
  function/ROL_LinearOperator.hpp
  function/ROL_ProjectedHessian.hpp
  function/ROL_ProjectedPreconditioner.hpp
  function/ROL_PrimalDualHessian.hpp
  function/ROL_PrimalDualPreconditioner.hpp
)

#
# step headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/step)
APPEND_SET(HEADERS
  step/ROL_Step.hpp
  step/ROL_PrimalDualActiveSetStep.hpp
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
  step/ROL_NonlinearCG.hpp
  step/ROL_CompositeStepSQP.hpp
  step/ROL_IterationScaling.hpp
  step/ROL_BackTracking.hpp
  step/ROL_CubicInterp.hpp
  step/ROL_Bisection.hpp
  step/ROL_GoldenSection.hpp
  step/ROL_Brents.hpp
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
# sol headers
#
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
INCLUDE_DIRECTORIES(${DIR}/function
                    ${DIR}/vector
                    ${DIR}/sampler
                    ${DIR}/sampler/sparse_grids
                   )
APPEND_GLOB(HEADERS ${DIR}/function/*.hpp 
                    ${DIR}/vector/*.hpp
                    ${DIR}/sampler/*.hpp
                    ${DIR}/sampler/sparse_grids/*.hpp
                   )

#
# zoo headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/zoo)
APPEND_SET(HEADERS
  zoo/ROL_Types.hpp
  zoo/ROL_HelperFunctions.hpp
  zoo/ROL_TestObjectives.hpp
  zoo/ROL_Rosenbrock.hpp
  zoo/ROL_FreudensteinRoth.hpp
  zoo/ROL_Beale.hpp
  zoo/ROL_Powell.hpp
  zoo/ROL_SumOfSquares.hpp
  zoo/ROL_LeastSquares.hpp
  zoo/ROL_PoissonControl.hpp
  zoo/ROL_PoissonInversion.hpp
  zoo/ROL_HS1.hpp
  zoo/ROL_HS2.hpp
  zoo/ROL_HS3.hpp
  zoo/ROL_HS4.hpp
  zoo/ROL_HS5.hpp
  zoo/ROL_HS25.hpp
  zoo/ROL_HS38.hpp
  zoo/ROL_HS45.hpp
  zoo/ROL_BVP.hpp
  zoo/ROL_Zakharov.hpp
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
