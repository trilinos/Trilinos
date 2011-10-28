INCLUDE(PackageLibraryMacros)

#
# Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# Define the header and source files (and directories)
#

SET(DIR  ${CMAKE_CURRENT_SOURCE_DIR})
SET(BDIR ${CMAKE_CURRENT_BINARY_DIR})

INCLUDE_DIRECTORIES(${DIR})
INCLUDE_DIRECTORIES(${DIR}/Headers)
INCLUDE_DIRECTORIES(${DIR}/MueCentral)
INCLUDE_DIRECTORIES(${DIR}/MueCentral/BaseClass)
INCLUDE_DIRECTORIES(${DIR}/MueCentral/Level)
INCLUDE_DIRECTORIES(${DIR}/Smoothers)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Smoothed-Aggregation)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Petrov-Galerkin-SA)
INCLUDE_DIRECTORIES(${DIR}/Transfers/Generic)
INCLUDE_DIRECTORIES(${DIR}/Transfers/BaseClass)
INCLUDE_DIRECTORIES(${DIR}/Graph)

SET(HEADERS "")
APPEND_GLOB(HEADERS ${BDIR}/${PACKAGE_NAME}_config.hpp)
APPEND_GLOB(HEADERS ${DIR}/Headers/${PACKAGE_NAME}_ConfigDefs.hpp)
APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(HEADERS */*/*.hpp)

SET(SOURCES "")
APPEND_GLOB(SOURCES *.cpp)
APPEND_GLOB(SOURCES */*.cpp)
APPEND_GLOB(SOURCES */*/*.cpp)

#
# Explicit instantiation
#

# Header file
PACKAGE_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/Graph)
INCLUDE_DIRECTORIES(${BDIR}/Graph)
APPEND_GLOB(HEADERS ${BDIR}/*.hpp)

# Cpp file
IF(${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION)
  CONFIGURE_FILE(${DIR}/../utils/ExplicitInstantiation/${PACKAGE_NAME}_ExplicitInstantiation.hpp.in ${BDIR}/MueLu_ExplicitInstantiation.hpp)
  APPEND_GLOB(SOURCES ${DIR}/../utils/ExplicitInstantiation/*.cpp)
ENDIF()

#
# Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
#MESSAGE("INCLUDE_DIRECTORIES = ")
#get_directory_property(output INCLUDE_DIRECTORIES)
#message(STATUS ${output})
