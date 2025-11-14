
INCLUDE(TribitsLibraryMacros)
INCLUDE(ROLUtils)

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

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

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

GET_PROPERTY( BLAS_STRING GLOBAL PROPERTY BLAS_IMPL )
GET_PROPERTY( BLAS_DIR    GLOBAL PROPERTY BLAS_DIR  )

GET_PROPERTY( MPI_STRING GLOBAL PROPERTY MPI_IMPL )
GET_PROPERTY( MPI_DIR    GLOBAL PROPERTY MPI_DIR  )

MESSAGE( "-- ROL has been configured to use the build options:"              )
MESSAGE( "--   ROL::Ptr           is implemented by ${PTR_STRING}"           )
MESSAGE( "--   ROL::ParameterList is implemented by ${PARAMETERLIST_STRING}" )
MESSAGE( "--   ROL::stacktrace    is implemented by ${STACKTRACE_STRING}"    )
MESSAGE( "--   ROL::LinearAlgebra is implemented by ${LA_STRING}"            )
MESSAGE( "--   ROL::LAPACK        is implemented by ${LAPACK_STRING}"        )
MESSAGE( "--   ROL::BLAS          is implemented by ${BLAS_STRING}"          )
MESSAGE( "--   ROL::MPI           is implemented by ${MPI_STRING}"           )

TRIBITS_INCLUDE_DIRECTORIES( ${PTR_DIR} )
APPEND_GLOB( HEADERS ${PTR_DIR}/*.hpp )

TRIBITS_INCLUDE_DIRECTORIES( ${PARAMETERLIST_DIR} )
APPEND_GLOB( HEADERS ${PARAMETERLIST_DIR}/*.hpp )

TRIBITS_INCLUDE_DIRECTORIES( ${STACKTRACE_DIR} )
APPEND_GLOB( HEADERS ${STACKTRACE_DIR}/*.hpp )

TRIBITS_INCLUDE_DIRECTORIES( ${LA_DIR} )
APPEND_GLOB( HEADERS ${LA_DIR}/*.hpp )

TRIBITS_INCLUDE_DIRECTORIES( ${LAPACK_DIR} )
APPEND_GLOB( HEADERS ${LAPACK_DIR}/*.hpp )

TRIBITS_INCLUDE_DIRECTORIES( ${BLAS_DIR} )
APPEND_GLOB( HEADERS ${BLAS_DIR}/*.hpp )

ROL_INCLUDE_DIRECTORIES( ${MPI_DIR} )
APPEND_GLOB( HEADERS ${MPI_DIR}/*.hpp )

#
# vector headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/vector)
TRIBITS_INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# function headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/function)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/dynamic
                    ${DIR}/boundconstraint
                    ${DIR}/constraint
                    ${DIR}/nlls
                    ${DIR}/polyproj
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
                    ${DIR}/polyproj/*.hpp
                    ${DIR}/objective/*.hpp
                    ${DIR}/operator/*.hpp
                    ${DIR}/simopt/*.hpp
                    ${DIR}/sketching/*.hpp
                    ${DIR}/std/*.hpp)

#
# step headers
#
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/step)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
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

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/algorithm)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/TypeB
                    ${DIR}/TypeE
                    ${DIR}/TypeG
                    ${DIR}/TypeU
										${DIR}/TypeP
                    ${DIR}/TypeB/pqn
                    ${DIR}/TypeG/augmentedlagrangian/
                    ${DIR}/TypeG/fletcher
                    ${DIR}/TypeG/interiorpoint
                    ${DIR}/TypeG/moreauyosida
                    ${DIR}/TypeG/stabilizedlcl
                    ${DIR}/TypeU/bundle
                    ${DIR}/TypeU/linesearch
                    ${DIR}/TypeU/linesearch/descent
                    ${DIR}/TypeU/trustregion)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/TypeB/*.hpp
                    ${DIR}/TypeE/*.hpp
                    ${DIR}/TypeG/*.hpp
                    ${DIR}/TypeU/*.hpp
										${DIR}/TypeP/*.hpp
                    ${DIR}/TypeB/pqn/*.hpp
                    ${DIR}/TypeG/augmentedlagrangian/*.hpp
                    ${DIR}/TypeG/fletcher/*.hpp
                    ${DIR}/TypeG/interiorpoint/*.hpp
                    ${DIR}/TypeG/moreauyosida/*.hpp
                    ${DIR}/TypeG/stabilizedlcl/*.hpp
                    ${DIR}/TypeU/bundle/*.hpp
                    ${DIR}/TypeU/linesearch/*.hpp
                    ${DIR}/TypeU/linesearch/descent/*.hpp
                    ${DIR}/TypeU/trustregion/*.hpp)

#
# parameterlist headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/parameterlist)
TRIBITS_INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# status test headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/status)
TRIBITS_INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# sol headers
#
TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/sol)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}/function
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
                    ${DIR}/algorithm/PrimalDual
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

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/elementwise)
TRIBITS_INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# utils headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/utils)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/function_bindings)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/function_bindings/*.hpp)

#
# zoo headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/zoo)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/testproblems)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/testproblems/*.hpp)
#
# oed headers
#

TRIBITS_SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/oed)
TRIBITS_INCLUDE_DIRECTORIES(${DIR}
                    ${DIR}/constraint
                    ${DIR}/factors
                    ${DIR}/objective
                    ${DIR}/objective/het
                    ${DIR}/objective/hom
                    ${DIR}/utilities)
APPEND_GLOB(HEADERS ${DIR}/*.hpp
                    ${DIR}/constraint/*.hpp
                    ${DIR}/factors/*.hpp
                    ${DIR}/objective/*.hpp
                    ${DIR}/objective/het/*.hpp
                    ${DIR}/objective/hom/*.hpp
                    ${DIR}/utilities/*.hpp)

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
