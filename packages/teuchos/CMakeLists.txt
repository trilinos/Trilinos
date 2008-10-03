# $Header$ 

PROJECT(Teuchos)

SET(PROJECT_DIR_NAME teuchos)
# 2008/08/08: ejguill: ToDO: Figure out PROJECT_DIR_NAME automatically from
# CMAKE_CURRENT_SOURCE_DIR

SET(${PROJECT_NAME}_VERSION "1.5d")

OPTION(${PROJECT_NAME}_ENABLE_COMPLEX
  "Enable Teuchos features for complex arithmetic."
  ON
  )

OPTION(${PROJECT_NAME}_ENABLE_DEBUG
  "Enable a host of runtime debug checking."
  OFF
  )

OPTION(${PROJECT_NAME}_ENABLE_GCC_DEMANGLE
  "Enable Teuchos calls to GCC demangle library function."
  OFF
  )

# 2008/10/02: rabartl: ToDo: Above: I need to add a macro like
# TRILINOS_DECLARE_OPTION(OPTION_NAME DEFINE_NAME DOCUMENTATION)

ADD_SUBDIRECTORY(src)

IF(${PROJECT_NAME}_ENABLE_TESTS)
  ADD_SUBDIRECTORY(test)
ENDIF()

IF(${PROJECT_NAME}_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(example)
ENDIF()

# Add targets runtests-serial-PROJECT and runtests-mpi-PROJECT
INCLUDE(TrilinosRuntestsTargets)
TRILINOS_ADD_CUSTOM_TEST_TARGETS()

LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_INCLUDE_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARY_DIRS)
LIST(REMOVE_DUPLICATES ${PROJECT_NAME}_LIBRARIES)

# 2008/09/15: rabartl: Above, we should write a function or macro to
# automatically clean up the duplicates?

##############################################################################
# Configure <PACKAGE>Config.cmake for the build tree.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.cmake.in
               ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake)
##############################################################################

# 2008/09/15: rabartl: Above, if this CONFIGURE_FILE(...) function is the same
# for every package why not just wrap it as a function that gets called?

# 2008/09/25: rabartl: Above, I don't think we even need a PackageConfig.cmake
# file
