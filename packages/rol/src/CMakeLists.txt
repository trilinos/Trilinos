
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

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/vector)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# function headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/function)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

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

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/elementwise)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

#
# zoo headers
#

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/zoo)
INCLUDE_DIRECTORIES(${DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)

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
