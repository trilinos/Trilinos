# $Header$ 

PROJECT(Teuchos)

SET(${PROJECT_NAME}_VERSION "1.5d")

# 2008/06/09: rabartl: Above: At the package level, I would prefer that you
# *not* use ${PROJECT_NAME} and instead just use Teuchos, or Epetra etc.
# Using the real package name makes greps more powerful and the concept of a
# PROJECT equating to a package is confusing.

# 2008/06/09: rabartl: I don't understand how a 'PROJECT' 'Trilinos' can be
# defined at the top Trilinos level and then redefined in each Trilinos
# package.

OPTION(${PROJECT_NAME}_ENABLE_TESTS
	"Build Teuchos tests."
	OFF
	)

OPTION(${PROJECT_NAME}_ENABLE_EXAMPLES
	"Build Teuchos examples."
	OFF
	)

OPTION(${PROJECT_NAME}_ENABLE_EXTENDED
	"Enable experimental Teuchos features."
	ON
	)

OPTION(${PROJECT_NAME}_ENABLE_COMPLEX
	"Enable Teuchos features for complex arithmetic."
	ON
	)

ADD_SUBDIRECTORY(src)

IF(${PROJECT_NAME}_ENABLE_TESTS)
	ADD_SUBDIRECTORY(test)
ENDIF(${PROJECT_NAME}_ENABLE_TESTS)

IF(${PROJECT_NAME}_ENABLE_EXAMPLES)
	ADD_SUBDIRECTORY(example)
ENDIF(${PROJECT_NAME}_ENABLE_EXAMPLES)

LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARIES)

##############################################################################
# Configure <PACKAGE>Config.cmake for the build tree.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
##############################################################################
