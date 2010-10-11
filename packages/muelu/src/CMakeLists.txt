
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.hpp)

#
# B) Define the header and source files (and directories)
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.hpp
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/MueCentral)

APPEND_SET(HEADERS
  MueCentral/MueLu_Version.hpp
  MueCentral/MueLu_ConfigDefs.hpp
  MueCentral/MueLu_Level.hpp
  Operator/MueLu_Operator.hpp
  ) 

APPEND_SET(SOURCES
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  muelu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

