
INCLUDE(TribitsLibraryMacros)
INCLUDE(TrilinosCreateClientTemplateHeaders)

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

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_GLOB(SOURCES ${DIR}/*.cpp)
TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR})

IF (${PACKAGE_NAME}_ENABLE_Experimental)
IF (${PACKAGE_NAME}_ENABLE_Boost)
IF (${PACKAGE_NAME}_ENABLE_Amesos2)
IF (${PACKAGE_NAME}_ENABLE_Cholmod)

TRILINOS_CREATE_CLIENT_TEMPLATE_HEADERS(${DIR}/supportgraph)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/supportgraph)

SET(HEADERS ${HEADERS}
supportgraph/Ifpack2_SupportGraph_decl.hpp
supportgraph/Ifpack2_SupportGraph_def.hpp
)
SET(SOURCES ${SOURCES}
supportgraph/Ifpack2_SupportGraph.cpp
)

ENDIF()
ENDIF()
ENDIF()
ENDIF()

# Must glob the binary dir last to get all of the auto-generated headers
SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_BINARY_DIR})
APPEND_GLOB(HEADERS ${DIR}/*.hpp)
APPEND_SET(HEADERS ${DIR}/${PACKAGE_NAME}_config.h )

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  ifpack2
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

