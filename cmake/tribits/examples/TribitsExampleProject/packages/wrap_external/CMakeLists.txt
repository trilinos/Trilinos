#
# Package: WrapExternal
#

# The purpose of this package is to demonstrate how to wrap an externally
# configured and built piece of software.  The external software once built
# provides libraries that need to be cast as TriBITS CMake libraries.  To
# make this even more interesting, the external software has dependencies on
# upstream TriBITS packages that it pulls in through a TriBITS export
# makefile.
#
# An overview of the basic tasks required are:
#
# a) Enumerate the upstream packages and TPLs needed (this is done in the
# cmake/Dependencies.cmake file).
#
# b) Generate an export makefile for the upstream packages and TPLs
#
# c) Call the configure script for the external package passing the compilers,
# compiler flags, and a pointer to the export makefile (which is just used to
# lift the include dirs and libs).
#
# d) Define library targets for the external project and define a build rule
# for generating those libraries given the list of header and source files in
# the external project.

TRIBITS_PACKAGE(WrapExternal)

#
# A) Determine if any downstream changes require a rebuild and deal with it
#

INCLUDE(TribitsFindMostRecentFileTimestamp)

ADVANCED_SET(${PACKAGE_NAME}_SHOW_MOST_RECENT_FILES   FALSE
  CACHE  BOOL
  "Set to TRUE to see most recent modified files.")
IF (${PACKAGE_NAME}_SHOW_MOST_RECENT_FILES)
  SET(SHOW_MOST_RECENT_FILES_ARG SHOW_MOST_RECENT_FILES)
ELSE()
  SET(SHOW_MOST_RECENT_FILES_ARG)
ENDIF()

TRIBITS_DETERMINE_IF_CURRENT_PACKAGE_NEEDS_REBUILT(
  ${SHOW_MOST_RECENT_FILES_ARG}
  SHOW_OVERALL_MOST_RECENT_FILES
  CURRENT_PACKAGE_OUT_OF_DATE_OUT ${PACKAGE_NAME}_BULID_IS_OUT_OF_DATE
  )

# If there is anything out of date, then just blow away the external build
IF (${PACKAGE_NAME}_BULID_IS_OUT_OF_DATE)
  MESSAGE(
    "\nBlowing away ${PACKAGE_NAME} build dir external_func/ so it will build from scratch!"
    )
  EXECUTE_PROCESS(
    COMMAND rm -rf "${PACKAGE_BINARY_DIR}/external_func/"
    )
ENDIF()
# NOTE: There may be other ways to address the out of date build other than to
# just blow the whole thing away but some external packages can be pretty
# complex and the only safe way to ensure that the rebuild is correct is to
# start from scratch.

#
# B) Write the export makefile that will be used by the external project
#

SET(EXPORT_MAKKEFILE "${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.TribExProj")
TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(
  PACKAGE_NAME ${PACKAGE_NAME}
  EXPORT_FILE_VAR_PREFIX TribitsExProj
  WRITE_EXPORT_MAKEFILE "${EXPORT_MAKKEFILE}"
  )

SET(EXTERNAL_FUNC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_func)
SET(EXTERNAL_FUNC_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/external_func)
SET(EXTERNAL_FUNC_LIB_FILE ${EXTERNAL_FUNC_BINARY_DIR}/libexternal_func.a)

FILE(MAKE_DIRECTORY ${EXTERNAL_FUNC_BINARY_DIR})

#
# C) Do configuration of the external project
#

EXECUTE_PROCESS(
  COMMAND ${PYTHON_EXECUTABLE} ${EXTERNAL_FUNC_SOURCE_DIR}/configure.py
    --with-export-makefile=${EXPORT_MAKKEFILE}
    --src-dir=${EXTERNAL_FUNC_SOURCE_DIR}
    --build-dir=${EXTERNAL_FUNC_BINARY_DIR}
  )

#
# D) Define a custom build rule and target to create external_func library
#

ADD_CUSTOM_COMMAND(
  OUTPUT ${EXTERNAL_FUNC_LIB_FILE}
  DEPENDS ${EXTERNAL_FUNC_SOURCE_DIR}/external_func.hpp
    ${EXTERNAL_FUNC_SOURCE_DIR}/external_func.cpp
  COMMAND make ${CTEST_BUILD_FLAGS}
  WORKING_DIRECTORY ${EXTERNAL_FUNC_BINARY_DIR}
  )
ADD_CUSTOM_TARGET( build_external_func
  DEPENDS ${EXTERNAL_FUNC_LIB_FILE} )
# NOTE: You have to create a custom target associated with the generated
# library.  You can't just use the custom command.
# NOTE: Above we list the source and header files that we know that if changed
# should trigger a rebuild of the library.  We would not likely know this with
# a more complex external project.

#
# E) Define the imported external library
#
# Below, does what TRIBITS_ADD_IMPORTED_LIBRARY() would do automatically!

# E.1) Create an imported library target and set up the dependencies to it
# custom target and command.
ADD_LIBRARY(external_func STATIC IMPORTED GLOBAL)
SET_PROPERTY(TARGET external_func PROPERTY IMPORTED_LOCATION
  ${EXTERNAL_FUNC_LIB_FILE})
ADD_DEPENDENCIES(external_func build_external_func)
# NOTE: Above, you have to use the custom target associated with the
# generation command and add it as a dependency of the imported library
# target.

# E.2) Make sure before we build the external library, we first build the
# libraries.
ADD_DEPENDENCIES(build_external_func pws_a)
# NOTE: You have to put the lib dependencies on build target, not the imported
# library target!

# E.3) Update the TriBITS variables
APPEND_SET(${PACKAGE_NAME}_LIB_TARGETS external_func)
GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES external_func pws_a)
GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${EXTERNAL_FUNC_SOURCE_DIR})
GLOBAL_SET(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES ON)
INCLUDE_DIRECTORIES(${EXTERNAL_FUNC_SOURCE_DIR})
# NOTE: Above, you have to add the upstream dependent libraries to the current
# package's list of libraries because you can't link to an importing lib with
# link_target_libraries() :-(

#
# F) Add an executable and test to show that it works!
#

IF (${PACKAGE_NAME}_ENABLE_TESTS)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(run_external_func
    SOURCES run_external_func.cpp
    PASS_REGULAR_EXPRESSION "external_func A no_deps"
    NUM_MPI_PROCS 1
    )

ENDIF()

# Run this if you need to ensure that the external library is build before you
# build this object file!
#SET_SOURCE_FILES_PROPERTIES(
#  run_external_func.cpp
#  PROPERTIES OBJECT_DEPENDS ${EXTERNAL_FUNC_LIB_FILE}
#  )

# Add a target to clean up the external function
ADD_CUSTOM_TARGET( clean_external_func
  COMMAND rm ${EXTERNAL_FUNC_LIB_FILE})
# NOTE: that we really only need to remove the final library.  The external build
# command will automatically know if the object file needs to be updated as
# well.

TRIBITS_PACKAGE_POSTPROCESS()

# LocalWords:  WrapExternal TriBITS CMake makefile TPLs
