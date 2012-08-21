
INCLUDE(TribitsLibraryMacros)


#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


SET(CURRENT_PACKAGE ShyLU)
CONFIGURE_FILE("${${PROJECT_NAME}_SOURCE_DIR}/packages/Trilinos_DLLExportMacro.h.in"
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_DLLExportMacro.h)


#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")


# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)
ASSERT_DEFINED(PROJECT_BINARY_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

#
# Include packages/src directory
#
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  shylu_debug_manager.hpp
  Ifpack_ShyLU.h
  Ifpack_Amesos_Schur.h
  shylu_util.h
  shylu_symbolic.h
  shylu.h
  shylu_probing_operator.h
  shylu_local_schur_operator.h
  shylu_internal.h
  )

# Set sources

APPEND_SET(SOURCES
  shylu_factor.cpp
  Ifpack_ShyLU.cpp
  Ifpack_Amesos_Schur.cpp
  shylu_probing_operator.cpp
  shylu_local_schur_operator.cpp
  shylu_util.cpp
  shylu_solve.cpp
  shylu_schur.cpp
  )

TRIBITS_ADD_LIBRARY(shylu
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
