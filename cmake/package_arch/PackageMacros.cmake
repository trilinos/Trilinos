INCLUDE(PackageSetupCompilerFlags)
INCLUDE(PackageWritePackageConfig)

INCLUDE(ParseVariableArguments)
INCLUDE(GlobalNullSet)
INCLUDE(AppendGlobalSet)
INCLUDE(PrintVar)
INCLUDE(PrependSet)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(AddOptionAndDefine)


#
# PACKAGE(...): Macro called at the very beginning of a ${PROJECT_NAME}
# package's top-level CMakeLists.txt file.
#
# Usage is:
#
#   PACKAGE(
#     <packageName>
#     [ENABLE_SHADOWING_WARNINGS]
#     [DISABLE_STRONG_WARNINGS]
#     [CLEANED]
#     [DISABLE_CIRCULAR_REF_DETECTION_FAILURE]
#     )
#
# The arguments are:
#
#   <packageName>
#
#     Gives the name of the Package, mostly just for checking and
#     documentation purposes.
#
#   ENABLE_SHADOWING_WARNINGS
#
#     If specified, then shadowing warnings will be turned on for supported
#     platforms/compilers.  The default is for shadowing warnings to be turned
#     off.  Note that this can be overridden globally by setting the cache
#     variable ${PROJECT_NAME}_ENABLE_SHADOWING_WARNIGNS.
#
#   DISABLE_STRONG_WARNINGS
#
#     If specified, then all strong warnings will be turned off, if they are
#     not already turned off by global cache varaibles.
#
#   CLEANED
#
#     If specified, then warnings will be promoted to errors for all defined
#     warnings.
#
#   DISABLE_CIRCULAR_REF_DETECTION_FAILURE
#
#     If specified, then the standard grep looking for RCPNode circular
#     references that causes tests to fail will be disabled.  Note that if
#     these warnings are being produced then it means that the test is leaking
#     memory and user like may also be leaking memory.
#

MACRO(PACKAGE PACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE: ${PACKAGE_NAME_IN}")
  ENDIF()
   
  #
  # A) Parse the input arguments
  #

  PARSE_ARGUMENTS(
    #prefix
    PARSE
    #lists
    ""
    #options
    "CLEANED;ENABLE_SHADOWING_WARNINGS;DISABLE_STRONG_WARNINGS;DISABLE_CIRCULAR_REF_DETECTION_FAILURE"
    ${ARGN}
    )

  #
  # B) Assert that the global and local package names are the same!
  #

  IF (DEFINED PACKAGE_NAME_GLOBAL)
    IF (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME_GLOBAL})
      MESSAGE(FATAL_ERROR "Error, the package-defined package name"
        " '${PACKAGE_NAME_IN}' is not the same as the package name"
        " defined at the global level '${PACKAGE_NAME_GLOBAL}'")
    ENDIF()
  ENDIF()

  #
  # C) Set up the CMake support for this ${PROJECT_NAME} package and define some
  # top-level varaibles.
  #

  SET(PACKAGE_NAME ${PACKAGE_NAME_IN})
  STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UC)
  IF (${PACKAGE_NAME}_ENABLE_TESTS AND ${PACKAGE_NAME}_ENABLE_EXAMPLES)
    SET(EXTRA_STATUS_MSG " (tests, examples)") 
  ELSEIF (${PACKAGE_NAME}_ENABLE_TESTS)
    SET(EXTRA_STATUS_MSG " (tests)") 
  ELSEIF (${PACKAGE_NAME}_ENABLE_EXAMPLES)
    SET(EXTRA_STATUS_MSG " (examples)") 
  ELSE()
    SET(EXTRA_STATUS_MSG "") 
  ENDIF() 
  MESSAGE("Processing enabled package: ${PACKAGE_NAME}${EXTRA_STATUS_MSG}")

  # Write PACKAGE versions of common variables
  SET(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this ${PROJECT_NAME} package is in
  FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  STRING(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

  # Set up the compile flags for the package
  PACKAGE_SETUP_COMPILER_FLASGS()

  # Set up circular reference detection test failure
  IF (PARSE_DISABLE_CIRCULAR_REF_DETECTION_FAILURE)
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF)
  ELSE()
    SET(${PACKAGE_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE ON)
  ENDIF()

  #
  # D) Define package linkage varaibles
  #

  GLOBAL_NULL_SET(${PACKAGE_NAME}_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_LIBRARIES)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_TEST_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_TEST_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_TEST_LIBRARIES)

  GLOBAL_NULL_SET(${PACKAGE_NAME}_LIB_TARGETS)
  GLOBAL_NULL_SET(${PACKAGE_NAME}_ALL_TARGETS)

ENDMACRO()


#
# Macro called to add a set of test directories for a package
#
# This macro only needs to be called from the top most CMakeList.txt file for
# which all subdirectories area all "tests".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# This macro defines hooks for inserting certain types of behavior in a
# uniform way.
#

MACRO(PACKAGE_ADD_TEST_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_TESTS)
    FOREACH(TEST_DIR ${ARGN})
      ADD_SUBDIRECTORY(${TEST_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Common options to add to a package
#


MACRO(PACKAGE_ADD_DEBUG_OPTION)
  ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_DEBUG
    HAVE_${PACKAGE_NAME_UC}_DEBUG
    "Enable a host of runtime debug checking."
    ${${PROJECT_NAME}_ENABLE_DEBUG}
    )
ENDMACRO()


MACRO(PACKAGE_ADD_SHOW_DEPRECATED_WARNINGS_OPTION)
  OPTION(
    ${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS
     "Show warnings about deprecated code in ${PACKAGE_NAME}"
    ${${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS}
    )
ENDMACRO()


MACRO(PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION)
  ADD_OPTION_AND_DEFINE(
    ${PACKAGE_NAME}_ENABLE_EXPLICIT_INSTANTIATION
    HAVE_${PACKAGE_NAME_UC}_EXPLICIT_INSTANTIATION
    "Enable the use of explicit template instantiation."
    ${${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION}
    )
ENDMACRO()


#
# Macro called to add a set of example directories for a package
#
# This macro only needs to be called from the top most CMakeList.txt file for
# which all subdirectories area all "examples".
#
# This macro can be called several times within a package and it will have the
# right effect.
#
# This macro defines hooks for inserting certain types of behavior in a
# uniform way.
#

MACRO(PACKAGE_ADD_EXAMPLE_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_EXAMPLES)
    FOREACH(EXAMPLE_DIR ${ARGN})
      ADD_SUBDIRECTORY(${EXAMPLE_DIR})
    ENDFOREACH()
  ENDIF()

ENDMACRO()


#
# Macro called at the very end of a ${PROJECT_NAME} package's top-level
# CMakeLists.txt file
#

MACRO(PACKAGE_POSTPROCESS)

  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_libs DEPENDS ${${PACKAGE_NAME}_LIB_TARGETS})
  ADD_CUSTOM_TARGET(${PACKAGE_NAME}_all DEPENDS ${${PACKAGE_NAME}_ALL_TARGETS})

  # Create the configure file so external projects can find packages with a
  # call to find_package(<package_name>)
  # This also creates the Makefile.export.* files.
  PACKAGE_WRITE_PACKAGE_CONFIG_FILE(${PACKAGE_NAME})

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  SET(${PACKAGE_NAME}_FINISHED_FIRST_CONFIGURE TRUE
    CACHE INTERNAL "")

ENDMACRO()


#
# Append the local package's cmake directory in order to help pull in 
# configure-time testing macros
#

PREPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
