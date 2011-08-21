
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "MueCentral/MueLu_Level.cpp")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Headers/MueLu_ConfigDefs.hpp
  )

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Headers)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/MueCentral)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Smoothers)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Transfers/Smoothed-Aggregation)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Transfers/Generic)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Transfers/BaseClass)
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/Graph)


APPEND_GLOB(HEADERS */*.hpp)
APPEND_GLOB(HEADERS */*/*.hpp)

# for debugging
#MESSAGE("")
#MESSAGE("**************** in muelu src ***********")
#MESSAGE("")
#MESSAGE("HEADERS = ${HEADERS}")
#MESSAGE("")
#MESSAGE("INCLUDE_DIRECTORIES = ")
#get_directory_property(output INCLUDE_DIRECTORIES)
#message(STATUS ${output})

APPEND_GLOB(SOURCES *.cpp)
APPEND_GLOB(SOURCES */*.cpp)
APPEND_GLOB(SOURCES */*/*.cpp)

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

