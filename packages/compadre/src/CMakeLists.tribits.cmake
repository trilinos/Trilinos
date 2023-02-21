#
# A) Package-specific configuration options
#

bob_config_header(${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h ${PACKAGE_NAME})

#
# B) Define the header and source files (and directories)
#

set(HEADERS "")
set(SOURCES "")

tribits_include_directories(${CMAKE_CURRENT_BINARY_DIR})

set(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_Config.h
  )

tribits_include_directories(${CMAKE_CURRENT_SOURCE_DIR})

append_glob(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
list(REMOVE_ITEM HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/Compadre_Manifold_Functions.hpp
  )
append_glob(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/basis/*.hpp)
append_glob(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/constraints/*.hpp)
append_glob(HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/tpl/*.hpp)
append_glob(SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

#
# C) Define the targets for package's library(s)
#

tribits_add_library(
  compadre
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# allows us to use flat directory includes when building, since that will be the file structure once installed
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/basis>)
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/constraints>)
target_include_directories(compadre PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/tpl>)
