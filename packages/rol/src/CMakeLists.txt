
INCLUDE(TribitsLibraryMacros)

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

#
# All headers
#

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
)

#
# compatibility headers (TODO: Move to compatibility/CMakeLists.txt)
#

GET_PROPERTY( PTR_STRING           GLOBAL PROPERTY PTR_IMPL )
GET_PROPERTY( PTR_DIR              GLOBAL PROPERTY PTR_DIR  )

GET_PROPERTY( STACKTRACE_STRING    GLOBAL PROPERTY STACKTRACE_IMPL )
GET_PROPERTY( STACKTRACE_DIR       GLOBAL PROPERTY STACKTRACE_DIR  )

GET_PROPERTY( PARAMETERLIST_STRING GLOBAL PROPERTY PARAMETERLIST_IMPL )
GET_PROPERTY( PARAMETERLIST_DIR    GLOBAL PROPERTY PARAMETERLIST_DIR  )

GET_PROPERTY( LA_STRING GLOBAL PROPERTY LA_IMPL )
GET_PROPERTY( LA_DIR    GLOBAL PROPERTY LA_DIR  )

GET_PROPERTY( LAPACK_STRING GLOBAL PROPERTY LAPACK_IMPL )
GET_PROPERTY( LAPACK_DIR    GLOBAL PROPERTY LAPACK_DIR  )

MESSAGE( "-- ROL has been configured to use the build options:"              )
MESSAGE( "--   ROL::Ptr           is implemented by ${PTR_STRING}"           )
MESSAGE( "--   ROL::ParameterList is implemented by ${PARAMETERLIST_STRING}" )
MESSAGE( "--   ROL::stacktrace    is implemented by ${STACKTRACE_STRING}"    )
MESSAGE( "--   ROL::LinearAlgebra is implemented by ${LA_STRING}"            )
MESSAGE( "--   ROL::LAPACK        is implemented by ${LAPACK_STRING}"        )

INCLUDE_DIRECTORIES( ${PTR_DIR} )
APPEND_GLOB( HEADERS ${PTR_DIR}/*.hpp )

INCLUDE_DIRECTORIES( ${PARAMETERLIST_DIR} )
APPEND_GLOB( HEADERS ${PARAMETERLIST_DIR}/*.hpp )

INCLUDE_DIRECTORIES( ${STACKTRACE_DIR} )
APPEND_GLOB( HEADERS ${STACKTRACE_DIR}/*.hpp )

INCLUDE_DIRECTORIES( ${LA_DIR} )
APPEND_GLOB( HEADERS ${LA_DIR}/*.hpp )

INCLUDE_DIRECTORIES( ${LAPACK_DIR} )
APPEND_GLOB( HEADERS ${LAPACK_DIR}/*.hpp )

#
# vector headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/vector)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# function headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/function)
INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/dynamic
                    ${DIR}/boundconstraint
                    ${DIR}/constraint
                    ${DIR}/nlls
                    ${DIR}/objective
                    ${DIR}/operator
                    ${DIR}/simopt
                    ${DIR}/sketching
                    ${DIR}/std)

APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/dynamic/*.hpp
                    ${DIR}/boundconstraint/*.hpp
                    ${DIR}/constraint/*.hpp
                    ${DIR}/nlls/*.hpp
                    ${DIR}/objective/*.hpp
                    ${DIR}/operator/*.hpp
                    ${DIR}/simopt/*.hpp
                    ${DIR}/sketching/*.hpp
                    ${DIR}/std/*.hpp)

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
                    ${DIR}/nonlinearcg
                    ${DIR}/fletcher
                   )
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/linesearch/*.hpp
                    ${DIR}/trustregion/*.hpp
                    ${DIR}/secant/*.hpp
                    ${DIR}/krylov/*.hpp
                    ${DIR}/bundle/*.hpp
                    ${DIR}/augmentedlagrangian/*.hpp
                    ${DIR}/moreauyosidapenalty/*.hpp
                    ${DIR}/interiorpoint/*.hpp
                    ${DIR}/nonlinearcg/*.hpp
		    ${DIR}/fletcher/*.hpp
                   )

#
# algorithm headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/algorithm)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# status test headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/status)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# sol headers
#
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
INCLUDE_DIRECTORIES(${DIR}/function
                    ${DIR}/function/distribution
                    ${DIR}/function/expectationquad
                    ${DIR}/function/randvarfunctional
                    ${DIR}/function/randvarfunctional/risk
                    ${DIR}/function/randvarfunctional/risk/fdivergence
                    ${DIR}/function/randvarfunctional/risk/spectral
                    ${DIR}/function/randvarfunctional/error
                    ${DIR}/function/randvarfunctional/regret
                    ${DIR}/function/randvarfunctional/deviation
                    ${DIR}/function/randvarfunctional/probability
                    ${DIR}/function/progressivehedging
#                    ${DIR}/function/riskmeasure
#                    ${DIR}/function/riskmeasure/quadrangle
#                    ${DIR}/function/riskmeasure/fdivergence
#                    ${DIR}/function/errormeasure
                    ${DIR}/algorithm
                    ${DIR}/vector
                    ${DIR}/sampler
                    ${DIR}/sampler/SROM
                    ${DIR}/status
                   )
APPEND_GLOB(HEADERS ${DIR}/function/*.hpp
                    ${DIR}/function/distribution/*.hpp
                    ${DIR}/function/riskmeasure/*.hpp
                    ${DIR}/function/expectationquad/*.hpp
                    ${DIR}/function/randvarfunctional/*.hpp
                    ${DIR}/function/randvarfunctional/risk/*.hpp
                    ${DIR}/function/randvarfunctional/risk/fdivergence/*.hpp
                    ${DIR}/function/randvarfunctional/risk/spectral/*.hpp
                    ${DIR}/function/randvarfunctional/error/*.hpp
                    ${DIR}/function/randvarfunctional/regret/*.hpp
                    ${DIR}/function/randvarfunctional/deviation/*.hpp
                    ${DIR}/function/randvarfunctional/probability/*.hpp
                    ${DIR}/function/progressivehedging/*.hpp
#                    ${DIR}/function/riskmeasure/quadrangle/*.hpp
#                    ${DIR}/function/riskmeasure/fdivergence/*.hpp
#                    ${DIR}/function/errormeasure/*.hpp
                    ${DIR}/algorithm/*.hpp
                    ${DIR}/vector/*.hpp
                    ${DIR}/sampler/*.hpp
                    ${DIR}/sampler/SROM/*.hpp
                    ${DIR}/status/*.hpp
                   )

#
# elementwise function headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/elementwise)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# utils headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils)
INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/function_bindings)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/function_bindings/*.hpp)

#
# zoo headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/zoo)
INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/testproblems)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/testproblems/*.hpp)

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
