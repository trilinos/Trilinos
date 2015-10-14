
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
  vector/ROL_PartitionedVector.hpp
  vector/ROL_ScaledStdVector.hpp
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
  function/ROL_LogBarrierObjective.hpp
)

#
# step headers
#
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/step)
INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/linesearch
                    ${DIR}/trustregion
                    ${DIR}/secant
                    ${DIR}/krylov
                    ${DIR}/bundle
                    ${DIR}/augmentedlagrangian
                    ${DIR}/moreauyosidapenalty
                    ${DIR}/interiorpoint
                   )
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/linesearch/*.hpp 
                    ${DIR}/trustregion/*.hpp
                    ${DIR}/secant/*.hpp
                    ${DIR}/krylov/*.hpp
                    ${DIR}/bundle/*.hpp
                    ${DIR}/augmentedlagrangian/*.hpp
                    ${DIR}/moreauyosidapenalty/*.hpp
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
  status/ROL_ConstraintStatusTest.hpp
  status/ROL_BundleStatusTest.hpp
  status/ROL_StatusTestFactory.hpp
  )

#
# sol headers
#
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
INCLUDE_DIRECTORIES(${DIR}/function
                    ${DIR}/function/distribution
                    ${DIR}/vector
                    ${DIR}/sampler
                    ${DIR}/sampler/SROM
                   )
APPEND_GLOB(HEADERS ${DIR}/function/*.hpp 
                    ${DIR}/function/distribution/*.hpp
                    ${DIR}/vector/*.hpp
                    ${DIR}/sampler/*.hpp
                    ${DIR}/sampler/SROM/*.hpp
                   )

#
# elementwise function headers
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/elementwise)
APPEND_SET(HEADERS
  elementwise/ROL_Elementwise_Function.hpp
  elementwise/ROL_Elementwise_Reduce.hpp
  elementwise/ROL_UnaryFunctions.hpp
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
  zoo/ROL_HS29.hpp
  zoo/ROL_HS32.hpp
  zoo/ROL_HS38.hpp
  zoo/ROL_HS45.hpp
  zoo/ROL_BVP.hpp
  zoo/ROL_Zakharov.hpp
  zoo/ROL_Minimax1.hpp
  zoo/ROL_Minimax2.hpp
  zoo/ROL_Minimax3.hpp
  zoo/ROL_ParameterListConverters.hpp
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
