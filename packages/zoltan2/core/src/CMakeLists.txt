#
# A) Package-specific configuration options
#

IF(NOT ${PROJECT_NAME}_ENABLE_OpenMP)
  SET(${PROJECT_NAME}_ENABLE_OpenMP OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE( Zoltan2_ENABLE_OpenMP
  HAVE_ZOLTAN2_OMP
  "Enable OpenMP support in Zoltan2."
  ${${PROJECT_NAME}_ENABLE_OpenMP}
)


TRIBITS_CONFIGURE_FILE(Zoltan2_config.h)

# Generate the header file that contains the Zoltan2 parameters.
# It is used only in compilation.  This header and the data that
# generates it need not be installed.

##KDD For now, disabled automatic generation of parameter header.
##KDD The generation didn't work when shared libraries were enabled.
##KDD Also, it would not work in a cross compiling environment.
##KDD The problems are due to limitations in CMake and TRIBITS.
##KDD After the Trilinos v11 release, we can revisit this issue.

##KDD SET(PARAM_HEADER ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_XML_Parameters.hpp)
##KDD SET(PARAM_DATA ${PACKAGE_SOURCE_DIR}/data/parameters.xml)
##KDD SET(PARAM_BUILDER_NAME xmlToHeaderDefinition)
##KDD SET(PARAM_BUILDER_SRC ${PACKAGE_SOURCE_DIR}/util/${PARAM_BUILDER_NAME}.cpp)
##KDD SET(PARAM_BUILDER_BIN ${CMAKE_CURRENT_BINARY_DIR}/${PARAM_BUILDER_NAME})

##KDD ADD_EXECUTABLE(${PARAM_BUILDER_NAME} ${PARAM_BUILDER_SRC})

##KDD ADD_CUSTOM_COMMAND(
##KDD   OUTPUT ${PARAM_HEADER}
##KDD   DEPENDS ${PARAM_BUILDER_BIN} ${PARAM_DATA}
##KDD   COMMAND ${PARAM_BUILDER_BIN}
##KDD   ARGS ${PARAM_DATA} ${PARAM_HEADER}
##KDD )

#
# B) Define the header files, the source files, and the include directories
#

SET(HEADERS "")
SET(SOURCES "")

##KDD APPEND_GLOB(HEADERS ${PARAM_HEADER})

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})

APPEND_GLOB(HEADERS ${DIR}/algorithms/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/partition/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/partition/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/order/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/order/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/color/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/color/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/match/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/match/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/map/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/map/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/algorithms/zoltan/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/algorithms/zoltan/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/environment/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/environment/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/input/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/input/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/models/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/models/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/problems/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/problems/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/util/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/util/*.cpp)

APPEND_GLOB(HEADERS ${DIR}/directory/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/directory/*.cpp)

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/Zoltan2_config.h
  ${PARAM_HEADER}
  )

# for Trilinos_version.h
ASSERT_DEFINED(PROJECT_BINARY_DIR)
INCLUDE_DIRECTORIES(${PROJECT_BINARY_DIR})

# for Zoltan2_config.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# for zoltan2 header files
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/partition)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/order)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/color)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/match)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/map)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/algorithms/zoltan)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/environment)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/input)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/models)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/problems)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/util)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/directory)

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  zoltan2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
