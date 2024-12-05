#
# Package: WrapExternal
#

# The purpose of this package is to demonstrate how to wrap an externally
# configured and built piece of software.  The external software once built
# provides libraries that need to be cast as TriBITS CMake libraries.  To
# make this even more interesting, the external software has dependencies on
# upstream TriBITS packages.

tribits_package(WrapExternal)

#
# A) Determine if any downstream changes require a rebuild and deal with it
#

include(TribitsFindMostRecentFileTimestamp)

advanced_set(${PACKAGE_NAME}_SHOW_MOST_RECENT_FILES   FALSE
  CACHE  BOOL
  "Set to TRUE to see most recent modified files.")
if (${PACKAGE_NAME}_SHOW_MOST_RECENT_FILES)
  set(SHOW_MOST_RECENT_FILES_ARG SHOW_MOST_RECENT_FILES)
else()
  set(SHOW_MOST_RECENT_FILES_ARG)
endif()

tribits_determine_if_current_package_needs_rebuilt(
  ${SHOW_MOST_RECENT_FILES_ARG}
  SHOW_OVERALL_MOST_RECENT_FILES
  CURRENT_PACKAGE_OUT_OF_DATE_OUT ${PACKAGE_NAME}_BULID_IS_OUT_OF_DATE
  )

# If there is anything out of date, then just blow away the external build
if (${PACKAGE_NAME}_BULID_IS_OUT_OF_DATE)
  message(
    "\nBlowing away ${PACKAGE_NAME} build dir external_func/ so it will build from scratch!"
    )
  execute_process(
    COMMAND rm -rf "${PACKAGE_BINARY_DIR}/external_func/"
    )
endif()
# NOTE: There may be other ways to address the out of date build other than to
# just blow the whole thing away but some external packages can be pretty
# complex and the only safe way to ensure that the rebuild is correct is to
# start from scratch.

#
# B) Get the arguments for the configure of external_func
#

string(TOUPPER "${CMAKE_BUILD_TYPE}" upperBuildType)

get_target_property(simpleCxx_IncludeDirs
  SimpleCxx::simplecxx  INTERFACE_INCLUDE_DIRECTORIES)

get_target_property(withSubpackagesA_IncludeDirs
  WithSubpackagesA::pws_a  INTERFACE_INCLUDE_DIRECTORIES)

set(includeDirsList ${withSubpackagesA_IncludeDirs} ${simpleCxx_IncludeDirs})

if (${PACKAGE_NAME}_ENABLE_MixedLang)
  get_target_property(mixedLang_IncludeDirs
    MixedLang::mixedlang  INTERFACE_INCLUDE_DIRECTORIES)
  list(PREPEND includeDirsList ${mixedLang_IncludeDirs})
endif()
list(JOIN includeDirsList " -I" includeDirsCompileOptions)
set(includeDirsCompileOptions "-I${includeDirsCompileOptions}")

#
# C) Do configuration of the external project as a build target
#

set(EXTERNAL_FUNC_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_func)
set(EXTERNAL_FUNC_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/external_func)
set(EXTERNAL_FUNC_LIB_FILE ${EXTERNAL_FUNC_BINARY_DIR}/libexternal_func.a)
set(EXTERNAL_FUNC_MAKEFILE ${EXTERNAL_FUNC_BINARY_DIR}/Makefile)

file(MAKE_DIRECTORY ${EXTERNAL_FUNC_BINARY_DIR})

add_custom_command(
  OUTPUT ${EXTERNAL_FUNC_MAKEFILE}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt
    ${EXTERNAL_FUNC_SOURCE_DIR}/configure.py
  COMMAND ${Python3_EXECUTABLE} ${EXTERNAL_FUNC_SOURCE_DIR}/configure.py
    "--cxx=${CMAKE_CXX_COMPILER}"
    "--cxx-flags=${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${upperBuildType}}"
    "--ar=${CMAKE_AR}"
    "--include-dirs=${includeDirsCompileOptions}"
    "--src-dir=${EXTERNAL_FUNC_SOURCE_DIR}"
    "--build-dir=${EXTERNAL_FUNC_BINARY_DIR}"
  )

# NOTE: Above, we have to run the configure.py script at build time after
# generation because ${includeDirsCompileOptions} contains generation
# expressions that are evaluated at generation time by the
# add_custom_command() call.

#
# D) Define a custom build rule and target to create external_func library
#

add_custom_command(
  OUTPUT ${EXTERNAL_FUNC_LIB_FILE}
  DEPENDS ${EXTERNAL_FUNC_SOURCE_DIR}/external_func.hpp
    ${EXTERNAL_FUNC_SOURCE_DIR}/external_func.cpp
    ${EXTERNAL_FUNC_MAKEFILE}
  COMMAND make ${CTEST_BUILD_FLAGS}
  WORKING_DIRECTORY ${EXTERNAL_FUNC_BINARY_DIR}
  )
add_custom_target( build_external_func
  DEPENDS ${EXTERNAL_FUNC_LIB_FILE} )
# NOTE: You have to create a custom target associated with the generated
# library.  You can't just use the custom command.
# NOTE: Above we list the source and header files that we know that if changed
# should trigger a rebuild of the library.  We would not likely know this with
# a more complex external project.

#
# E) Define the imported external library
#
# Below, does what tribits_add_imported_library() would do automatically!

# E.1) Create an imported library target and set up the dependencies to it
# custom target and command.
add_library(external_func STATIC IMPORTED GLOBAL)
set_property(TARGET external_func PROPERTY IMPORTED_LOCATION
  ${EXTERNAL_FUNC_LIB_FILE})
add_dependencies(external_func build_external_func)
# NOTE: Above, you have to use the custom target associated with the
# generation command and add it as a dependency of the imported library
# target.

# E.2) Make sure before we build the external library, we first build the
# libraries.
add_dependencies(build_external_func WithSubpackagesA::pws_a)
# NOTE: You have to put the lib dependencies on build target, not the imported
# library target!

# E.3) Update the TriBITS variables
append_set(${PACKAGE_NAME}_LIB_TARGETS external_func)
global_set(${PACKAGE_NAME}_LIBRARIES external_func WithSubpackagesA::pws_a)
global_set(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES ON)
tribits_include_directories(${EXTERNAL_FUNC_SOURCE_DIR})
# NOTE: Above, you have to add the upstream dependent libraries to the current
# package's list of libraries because you can't link to an importing lib with
# link_target_libraries() :-(

#
# F) Add an executable and test to show that it works!
#

if (${PACKAGE_NAME}_ENABLE_TESTS)

  tribits_add_executable_and_test(run_external_func
    SOURCES run_external_func.cpp
    PASS_REGULAR_EXPRESSION "external_func A ${EXPECTED_SIMPLECXX_AND_DEPS}"
    NUM_MPI_PROCS 1
    )

endif()

# Run this if you need to ensure that the external library is build before you
# build this object file!
#set_source_files_properties(
#  run_external_func.cpp
#  PROPERTIES OBJECT_DEPENDS ${EXTERNAL_FUNC_LIB_FILE}
#  )

# Add a target to clean up the external function
add_custom_target( clean_external_func
  COMMAND rm ${EXTERNAL_FUNC_LIB_FILE})
# NOTE: that we really only need to remove the final library.  The external build
# command will automatically know if the object file needs to be updated as
# well.

tribits_package_postprocess()

# LocalWords:  WrapExternal TriBITS CMake makefile TPLs
