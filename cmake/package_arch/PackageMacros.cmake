INCLUDE(ParseVariableArguments)
INCLUDE(GlobalNullSet)
INCLUDE(AppendGlobalSet)
INCLUDE(PrintVar)
INCLUDE(PrependSet)
INCLUDE(RemoveGlobalDuplicates)


#
# Macro called at the very beginning of a ${PROJECT_NAME} package's top-level
# CMakeLists.txt file
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
    "CLEANED"
    ${ARGN}
    )

  #
  # B) Assert that the global and local package names are the same!
  #

  IF (DEFINED PACKAGE_NAME_GLOBAL)
    IF (NOT ${PACKAGE_NAME_IN} STREQUAL ${PACKAGE_NAME_GLOBAL})
      MESSAGE(FATAL_ERROR "Error, the package-defined package name '${PACKAGE_NAME_IN}' is not the same as the package name defined at the global level '${PACKAGE_NAME_GLOBAL}'")
    ENDIF()
  ENDIF()

  #
  # C) Set up the CMake support for this ${PROJECT_NAME} package and define some
  # top-level varaibles.
  #

  SET(PACKAGE_NAME ${PACKAGE_NAME_IN})
  MESSAGE(STATUS "Processing enabled package: ${PACKAGE_NAME}")

  # Write PACKAGE versions of common variables
  SET(PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")
  SET(PACKAGE_BINARY_DIR "${CMAKE_CURRENT_BINARY_DIR}")

  # Get the name of the directory this ${PROJECT_NAME} package is in
  FILE(TO_CMAKE_PATH ${CMAKE_CURRENT_SOURCE_DIR} STANDARD_PACKAGE_SOURCE_DIR)
  STRING(REGEX REPLACE "/.+/(.+)" "\\1" PACKAGE_DIR_NAME "${STANDARD_PACKAGE_SOURCE_DIR}")

  # Set up for warnings as errors if requested

  ASSERT_DEFINED(PARSE_CLEANED)
	

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C ${PROJECT_NAME}_ENABLE_C_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS AND CMAKE_BUILD_TYPE)
    SET(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}
      " ${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS} ${CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE}}") 
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Setting up for C warnings as errors just in this package ...")
      PRINT_VAR(CMAKE_C_FLAGS_${CMAKE_BUILD_TYPE})
    ENDIF()
  ENDIF()

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX ${PROJECT_NAME}_ENABLE_CXX_DEBUG_COMPILE_FLAGS)
  IF (PARSE_CLEANED AND ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS AND CMAKE_BUILD_TYPE)
    SET(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}
      " ${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS} ${CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE}}") 
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Setting up for C++ warnings as errors just in this package ...")
      PRINT_VAR(CMAKE_CXX_FLAGS_${CMAKE_BUILD_TYPE})
    ENDIF()
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

  #
  # E) Define standard runtests targets for home-grown perl-based test harness
  #

  IF (${PROJECT_NAME}_ENABLE_NATIVE_TEST_HARNESS)
  
    ADD_CUSTOM_TARGET(
      ${PACKAGE_NAME}-runtests-serial
       ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
      --trilinos-dir=${TRILINOS_HOME_DIR}
      --comm=serial
      --build-dir=${TRILINOS_BUILD_DIR}
      --category=${TRILINOS_TEST_CATEGORY}
      --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
      --verbosity=1
      --packages=${PACKAGE_DIR_NAME}
      )

    IF (TPL_ENABLE_MPI)
    
      ADD_CUSTOM_TARGET(
        ${PACKAGE_NAME}-runtests-mpi
         ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
        --trilinos-dir=${TRILINOS_HOME_DIR}
        --comm=mpi
        --mpi-go="${TRILINOS_MPI_GO}"
        --max-proc=${MPIEXEC_MAX_NUMPROCS}
        --build-dir=${TRILINOS_BUILD_DIR}
        --category=${TRILINOS_TEST_CATEGORY}
        --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
        --verbosity=1
        --packages=${PACKAGE_DIR_NAME}
        )

    ENDIF()

  ENDIF()

  # 2008/12/23: rabartl: ToDo: Above: Get rid of these targets since
  # we will not need them very soon.

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
# Macro called to add a set of performance test directories for a package
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

MACRO(PACKAGE_ADD_PERFORMANCE_TEST_DIRECTORIES)

  IF(${PACKAGE_NAME}_ENABLE_PERFORMANCE_TESTS OR ${PACKAGE_NAME}_ENABLE_TESTS)
    FOREACH(TEST_DIR ${ARGN})
      ADD_SUBDIRECTORY(${TEST_DIR})
    ENDFOREACH()
  ENDIF()

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

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_POSTPROCESS: ${PACKAGE_NAME}")
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

ENDMACRO()


#
# Append the local package's cmake directory in order to help pull in 
# configure-time testing macros
#

PREPEND_SET(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
