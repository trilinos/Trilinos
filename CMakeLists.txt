

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)


#
# A) Basic top-level Trilinos stuff
#

MESSAGE("")
MESSAGE("Configuring Trilinos build directory")
MESSAGE("")

PROJECT(Trilinos NONE)

SET(Trilinos_VERSION "9.0d")
SET(TRILINOS_HOME_DIR ${CMAKE_CURRENT_SOURCE_DIR} CACHE INTERNAL "")
SET(TRILINOS_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR} CACHE INTERNAL "")
# Above, we put these in the cache so we can grep them out

SET( CMAKE_MODULE_PATH
   ${TRILINOS_HOME_DIR}/cmake
   ${TRILINOS_HOME_DIR}/cmake/utils
   ${TRILINOS_HOME_DIR}/cmake/package_arch
   ${TRILINOS_HOME_DIR}/cmake/config_tests
   )

INCLUDE(PackageArchGlobalMacros)
INCLUDE(TrilinosGlobalMacros)
INCLUDE(AdvancedSet)
INCLUDE(AdvancedOption)

OPTION(Trilinos_VERBOSE_CONFIGURE
  "Make the Trilinos configure process verbose." OFF)


#
# A.1) Run misc unit tests that don't need anything else
#

IF (Trilinos_INVOKE_TESTING_UNIT_TESTS)
  ADD_SUBDIRECTORY(cmake/TestingUnitTests/UnitTests)
  RETURN()
ENDIF()

#
# A.2) Set up other stuff
#

# Find Python executable which is needed for dependency file building
FIND_PACKAGE(PythonInterp)


#
# B) Define the list of all Trilinos packages and TPLs
#

# B.1) Packages

IF (NOT Trilinos_PACKAGES_FILE)
  SET(Trilinos_PACKAGES_FILE "TrilinosPackages")
ENDIF()
# NOTE: Above, I am not using a cache variable because I don't want
# the general user to even see this option!

MESSAGE("")
MESSAGE("Reading the list of packages from ${Trilinos_PACKAGES_FILE}.cmake ... ")
MESSAGE("")

INCLUDE(${Trilinos_PACKAGES_FILE})

PACKAGE_ARCH_PROCESS_PACKAGES_AND_DIRS_LISTS()

# B.2) TPLs

IF (NOT Trilinos_TPLS_FILE)
  SET(Trilinos_TPLS_FILE "TrilinosTPLs")
ENDIF()
# NOTE: Above, I am not using a cache variable because I don't want
# the general user to even see this option!

MESSAGE("")
MESSAGE("Reading the list of TPLs from ${Trilinos_TPLS_FILE}.cmake ... ")
MESSAGE("")

INCLUDE(${Trilinos_TPLS_FILE})

PACKAGE_ARCH_PROCESS_TPLS_LISTS()


#
# C) Process the package and TPL dependencies
#

PACKAGE_ARCH_READ_ALL_PACKAGE_DEPENDENCIES()


#
# D) Set up user options and global variables that will be used throughout
#

MESSAGE("")
MESSAGE("Setting up major user options ...")
MESSAGE("")

PACKAGE_ARCH_DEFINE_GLOBAL_OPTIONS()

IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
  SET(Trilinos_ENABLE_DEBUG_DEFAULT ON)
ELSE()
  SET(Trilinos_ENABLE_DEBUG_DEFAULT OFF)
ENDIF()
SET(Trilinos_ENABLE_DEBUG ${Trilinos_ENABLE_DEBUG_DEFAULT} CACHE BOOL
  "Enable debug checking for Trilinos packages.  Off by default unless CMAKE_BUILD_TYPE=\"DEBUG\"." )

ADVANCED_OPTION(Trilinos_ENABLE_EXPLICIT_INSTANTIATION
  "Enable explicit template instantiation in all package where it is supported."
  OFF )

ADVANCED_SET(Trilinos_DATA_DIR
  "NOTFOUND"
  CACHE PATH
  "Path TrilinosData directory to find more tests and other stuff" )

ADVANCED_OPTION(Trilinos_ENABLE_PERFORMANCE_TESTS
  "Enable performance tests in all packages that support them"
  OFF )

ADVANCED_OPTION(Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING
  "Shortcircut after dependency handling is complete"
  OFF )

IF (WIN32 AND NOT CYGWIN)
  SET(Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT "")
ELSE()
  SET(Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_HOME_DIR}/cmake/python/data/TrilinosPackageDependencies.xml" )
ENDIF()
ADVANCED_SET(Trilinos_DEPS_XML_OUTPUT_FILE
  "${Trilinos_DEPS_XML_OUTPUT_FILE_DEFAULT}"
  CACHE STRING
  "Output XML file containing Trilinos dependenices used by tools (if not empty)." )
# 2009/01/19: rabartl: Above: This file outputs just fine on MS
# Windows using MS Visual Studio but it causes the entire file to
# diff.  There must be something wrong with a newlines or something
# that is causing this.  If people are going to be doing real
# development work on MS Windows with MS Visual Studio, then we need
# to fix this so that the dependency files will get created and
# checked in correctly.  I will look into this later.

IF(Trilinos_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE)
  SET(Trilinos_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_HOME_DIR}/cmake/python/data/CDashSubprojectDependencies.xml" )
ELSE()
  SET(Trilinos_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT "")
ENDIF()
ADVANCED_SET(Trilinos_CDASH_DEPS_XML_OUTPUT_FILE
  "${Trilinos_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT}"
  CACHE STRING
  "Output XML file used by CDash in Trilinos-independent format (if not empty)." )

IF(Trilinos_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE)
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT
    "${TRILINOS_BUILD_DIR}/TrilinosPackageDependenciesTable.html" )
ELSE()
  SET(Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT "")
ENDIF()
ADVANCED_SET(Trilinos_DEPS_HTML_OUTPUT_FILE
  "${Trilinos_DEPS_HTML_OUTPUT_FILE_DEFAULT}"
  CACHE STRING
  "HTML Trilinos dependenices file that will be written to (if not empty)." )


#
# E) Apply logic to enable Trilinos packages and tests
#

PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
  "\nExplicitly enabled packages on input (by user)" ON FALSE)
PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
  "\nExplicitly disabled packages on input (by user or by default)" OFF FALSE)
PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
  "\nExplicitly enabled TPLs on input (by user)" ON FALSE)
PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
  "\nExplicitly disabled TPLs on input (by user or by default)" OFF FALSE)

PACKAGE_ARCH_ADJUST_PACKAGE_ENABLES()

PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
  "\nFinal set of enabled packages" ON FALSE)
PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
  "\nFinal set of non-enabled packages" OFF TRUE)
PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
  "\nFinal set of enabled TPLs" ON FALSE)
PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
  "\nFinal set of non-enabled TPLs" OFF TRUE)

IF (Trilinos_DEPS_XML_OUTPUT_FILE)
  IF (NOT IS_ABSOLUTE Trilinos_DEPS_XML_OUTPUT_FILE)
    SET(Trilinos_DEPS_XML_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${Trilinos_DEPS_XML_OUTPUT_FILE})
  ENDIF()
  MESSAGE("" )
  MESSAGE("Dumping the XML dependencies file ${Trilinos_DEPS_XML_OUTPUT_FILE} ..." )
  MESSAGE("")
  PACKAGE_ARCH_DUMP_DEPS_XML_FILE()
ENDIF()

IF (Trilinos_DEPS_HTML_OUTPUT_FILE)
  MESSAGE("" )
  MESSAGE("Dumping the HTML dependencies webpage file ${Trilinos_DEPS_HTML_OUTPUT_FILE} ..." )
  MESSAGE("" )
  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE}
      ${TRILINOS_HOME_DIR}/cmake/python/dump-package-dep-table.py
      --input-xml-deps-file=${Trilinos_DEPS_XML_OUTPUT_FILE}
      --output-html-deps-file=${Trilinos_DEPS_HTML_OUTPUT_FILE} )
ENDIF()

IF (Trilinos_CDASH_DEPS_XML_OUTPUT_FILE)
  MESSAGE("" )
  MESSAGE("Dumping the CDash XML dependencies file ${Trilinos_CDASH_DEPS_XML_OUTPUT_FILE} ..." )
  MESSAGE("" )
  EXECUTE_PROCESS(
    COMMAND ${PYTHON_EXECUTABLE}
      ${TRILINOS_HOME_DIR}/cmake/python/dump-cdash-deps-xml-file.py
      --input-xml-deps-file=${Trilinos_DEPS_XML_OUTPUT_FILE}
      --output-cdash-deps-xml-file=${Trilinos_CDASH_DEPS_XML_OUTPUT_FILE} )
ENDIF()

IF (Trilinos_SHORTCIRCUIT_AFTER_DEPENDENCY_HANDLING)
  MESSAGE("")
  MESSAGE(SEND_ERROR "Shortcircuiting after dependency tracking ...")
  RETURN()
ENDIF()


#
# F) Set up the environment on this computer
#

MESSAGE("")
MESSAGE("Probing the environment ...")
MESSAGE("")

PACKAGE_ARCH_PRE_SETUP_ENV()

INCLUDE(TrilinosProbeEnv)

PACKAGE_ARCH_POST_SETUP_ENV()


#
# G) Go get the information for all enabled TPLS
#

MESSAGE("")
MESSAGE("Getting information for all enabled TPLs ...")
MESSAGE("")

PACKAGE_ARCH_PROCESS_ENABLED_TPLS()


#
# H) Set up for testing with CTest and Trilinos test harness
#

MESSAGE("")
MESSAGE("Setting up testing support ...")
MESSAGE("")

INCLUDE(CTest)

TRILINOS_SETUP_TESTING_SUPPORT()


#
# I) Setup for packaging and distribution
#

TRILINOS_DEFINE_PACKAGING()


#
# J) Configure individual packages
# 

MESSAGE("")
MESSAGE("Configuring individual enabled Trilinos packages ...")
MESSAGE("")

PACKAGE_ARCH_CONFIGURE_ENABLED_PACKAGES()


#
# K) Install-related commands
#

CONFIGURE_FILE( ${CMAKE_CURRENT_SOURCE_DIR}/cmake/TrilinosConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/TrilinosConfig.cmake )
