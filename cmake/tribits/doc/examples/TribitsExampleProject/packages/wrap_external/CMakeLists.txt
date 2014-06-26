#
# Package: WrapExternal
#

# The purpose of this package is to demonstrate how to wrap an externally
# configured and built piece of software.  The external software once built
# provides libraries that need to be cast as TriBITS CMake librarires.  To
# make this even more interesting, the external software has dependancies on
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
# compiler flags, and a pointer to the export makefile (which they just use to
# grab the include dirs and libs out of).
#
# d) Define library targets for the external project and define a build rule
# for generating those libraries given the list of header and source files in
# the external project.

TRIBITS_PACKAGE(WrapExternal)

# A) Write the export makefile that will be used by the external project
SET(EXPORT_MAKKEFILE "${CMAKE_CURRENT_BINARY_DIR}/Makefile.export.TribExProj")
TRIBITS_WRITE_FLEXIBLE_PACKAGE_CLIENT_EXPORT_FILES(
  PACKAGE_NAME ${PACKAGE_NAME}
  EXPORT_FILE_VAR_PREFIX TribitsExProj
  WRITE_EXPORT_MAKLEFILE "${EXPORT_MAKKEFILE}"
  )

SET(EXTERNAL_FUNC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_func)
#PRINT_VAR(EXTERNAL_FUNC_SOURCE_DIR)
SET(EXTERNAL_FUNC_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/external_func)
SET(EXTERNAL_FUNC_LIB_FILE ${EXTERNAL_FUNC_BINARY_DIR}/libexternal_func.a)
#PRINT_VAR(EXTERNAL_FUNC_LIB_FILE)

FILE(MAKE_DIRECTORY ${EXTERNAL_FUNC_BINARY_DIR})

# B) Do configuration of the external project
EXECUTE_PROCESS(
  COMMAND ${PYTHON_EXECUTABLE} ${EXTERNAL_FUNC_SOURCE_DIR}/configure.py
    --with-export-makefile=${EXPORT_MAKKEFILE}
    --src-dir=${EXTERNAL_FUNC_SOURCE_DIR}
    --build-dir=${EXTERNAL_FUNC_BINARY_DIR}
  )

# C) Define a custom build rule and target to create exteranl_func library
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
  
# D) Add the imported library with TRIBITS_ADD_LIBRARY( ... IMPORTED ...)
#
# Below, I just manuall do what TRIBITS_ADD_LIBRARY() would do automatically.

# D.1) Create an imported library target and set up the depenancies
ADD_LIBRARY(exteranl_func STATIC IMPORTED GLBOAL) # GLOBAL
SET_PROPERTY(TARGET exteranl_func PROPERTY IMPORTED_LOCATION ${EXTERNAL_FUNC_LIB_FILE})
ADD_DEPENDENCIES(build_external_func pws_c)
ADD_DEPENDENCIES(exteranl_func build_external_func)
GLOBAL_SET(exteranl_func_IMPORTLIB_TARGET build_external_func)
# NOTE: Above, you have to use the custom target associated with the
# generation command and add it as a dependancy of the imported library
# target.

# D.2) Update the TriBITS varaibles
APPEND_SET(${PACKAGE_NAME}_LIB_TARGETS exteranl_func)
GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES exteranl_func pws_c)
GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${EXTERNAL_FUNC_SOURCE_DIR})
GLOBAL_SET(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES ON)
INCLUDE_DIRECTORIES(${EXTERNAL_FUNC_SOURCE_DIR})
# NOTE: Above, you have to add the upstream dependent libraries to the current
# package's list of libraries because you can't link to an importing lib with
# link_target_libraries() :-(

#
# E) Add an executable and test to show that it works!
#

TRIBITS_ADD_EXECUTABLE_AND_TEST(run_external_func
  SOURCES run_external_func.cpp
  DEPLIBS exteranl_func
  PASS_REGULAR_EXPRESSION "external_func C B A"
  NUM_MPI_PROCS 1
  )

ADD_CUSTOM_TARGET( clean_exteranl_func
  COMMAND rm external_func/external_func.o)

TRIBITS_PACKAGE_POSTPROCESS()
