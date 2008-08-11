# $Header$

# Project name
PROJECT(Epetra CXX C)

SET(PROJECT_DIR_NAME epetra)
# 2008/08/07: rabartl: ToDO: Figure out PROJECT_DIR_NAME automatically from
# CMAKE_CURRENT_SOURCE_DIR

# Setup package version
SET(${PROJECT_NAME}_VERSION "3.7d")

# Flag to enable tests
OPTION(${PROJECT_NAME}_ENABLE_TESTS
       "Build ${PROJECT_NAME} tests."
       OFF
)

# Flag to enable examples
OPTION(${PROJECT_NAME}_ENABLE_EXAMPLES
       "Build ${PROJECT_NAME} examples."
       OFF
)

# Flag to enable Fortran implementation code for performance
IF(WIN32 AND NOT CYGWIN)
	SET(${PROJECT_NAME}_ENABLE_FORTRAN_DEFAULT OFF)
ELSE(WIN32 AND NOT CYGWIN)
	SET(${PROJECT_NAME}_ENABLE_FORTRAN_DEFAULT ON)
ENDIF(WIN32 AND NOT CYGWIN)

OPTION(${PROJECT_NAME}_ENABLE_FORTRAN
	"Build Fortran library for this package."
	${${PROJECT_NAME}_ENABLE_FORTRAN_DEFAULT}
	)

SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR/cmake})

# Fortran related configuration
IF(${PROJECT_NAME}_ENABLE_FORTRAN)
  ENABLE_LANGUAGE(Fortran)
ENDIF(${PROJECT_NAME}_ENABLE_FORTRAN)

ADD_SUBDIRECTORY(src)

IF(${PROJECT_NAME}_ENABLE_TESTS)
  ADD_SUBDIRECTORY(test)
ENDIF(${PROJECT_NAME}_ENABLE_TESTS)

IF(${PROJECT_NAME}_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(example)
ENDIF(${PROJECT_NAME}_ENABLE_EXAMPLES)

INCLUDE(TrilinosRuntestsTargets)

LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARIES)

##############################################################################
# Configure <PACKAGE>Config.cmake for the build tree.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
##############################################################################

