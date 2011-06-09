
INCLUDE(PackageArchConstants)
INCLUDE(PackageArchProcessExtraExternalRepositoriesLists)
INCLUDE(PackageArchProcessPackagesAndDirsLists)
INCLUDE(PackageArchAdjustPackageEnables)
INCLUDE(PackageArchSetupMPI)
INCLUDE(PackageArchCategories)

INCLUDE(AddOptionAndDefine)
INCLUDE(AdvancedOption)
INCLUDE(AdvancedSet)
INCLUDE(AppendStringVar)
INCLUDE(AssertAndTouchDefined)
INCLUDE(CMakeBuildTypesList)
INCLUDE(FindListElement)
INCLUDE(GlobalNullSet)
INCLUDE(PrintNonemptyVar)
INCLUDE(PrintVar)
INCLUDE(Split)
INCLUDE(RemoveGlobalDuplicates)




#
# Define and option to include a file that reads in a bunch of options
#
#

MACRO(PACKAGE_ARCH_READ_IN_OPTIONS_FROM_FILE)


  SET( ${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE "" CACHE FILEPATH
    "Name of an optional file that is included first to define any cmake options with SET( ... CACHE ...) calls."
    )

  FOREACH (CONFIG_OPTS_FILE ${${PROJECT_NAME}_CONFIGURE_OPTIONS_FILE})
    MESSAGE("Reading in configuration options from ${CONFIG_OPTS_FILE} ...")
    INCLUDE(${CONFIG_OPTS_FILE})
  ENDFOREACH()


ENDMACRO()


#
# Define all of the standard global package architecture options.
#

MACRO(PACKAGE_ARCH_DEFINE_GLOBAL_OPTIONS)

  SET( ${PROJECT_NAME}_ENABLE_ALL_PACKAGES OFF CACHE BOOL
    "Enable all packages (Primary Stable and perhaps Secondary Stable packages)." )
  
  SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES ON CACHE BOOL
    "Recursively enable all optional packages for set of enabled packages." )
  
  ADVANCED_SET(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL
    "Recursively enable all packages that have required or optional dependencies for set of enabled packages." )
  
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_TESTS ""
    "Enable tests in all packages  (set to ON, OFF, or leave empty)." )
  
  SET_CACHE_ON_OFF_EMPTY(${PROJECT_NAME}_ENABLE_EXAMPLES ""
    "Enable examples in all packages  (set to ON, OFF, or leave empty).  If left empty, then this will be set to ON if ${PROJECT_NAME}_ENABLE_TESTS=ON" )
  
  IF (${PROJECT_NAME}_ENABLE_TESTS AND ${PROJECT_NAME}_ENABLE_EXAMPLES STREQUAL "")
    MESSAGE(STATUS "Setting ${PROJECT_NAME}_ENABLE_EXAMPLES=ON because ${PROJECT_NAME}_ENABLE_TESTS=ON")
    SET(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
    "Set to empty all package enables (set to OFF at end)." )
  
  ADVANCED_OPTION(${PROJECT_NAME}_ENABLE_C
    "Enable the C compiler and related code"
    ON )
  
  ADVANCED_OPTION(${PROJECT_NAME}_ENABLE_CXX
    "Enable the C++ compiler and related code"
    ON )
  
  IF(WIN32 AND NOT CYGWIN)
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT OFF)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_Fortran_DEFAULT ON)
  ENDIF()
  
  OPTION(${PROJECT_NAME}_ENABLE_Fortran
    "Enable the Fortran compiler and related code"
    ${${PROJECT_NAME}_ENABLE_Fortran_DEFAULT} )
  
  ADVANCED_SET(${PROJECT_NAME}_EXTRA_LINK_FLAGS ""
    CACHE STRING
    "Extra flags added to the end of every linked executable"
    )

  IF (CMAKE_BUILD_TYPE STREQUAL "DEBUG")
    SET(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT ON)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT OFF)
  ENDIF()
  SET(${PROJECT_NAME}_ENABLE_DEBUG ${${PROJECT_NAME}_ENABLE_DEBUG_DEFAULT} CACHE BOOL
    "Enable debug checking for ${PROJECT_NAME} packages.  Off by default unless CMAKE_BUILD_TYPE=\"DEBUG\"." )
  
  ADVANCED_SET(${PROJECT_NAME}_SHOW_DEPRECATED_WARNINGS ON
    CACHE BOOL
    "Show warnings about deprecated code"
    )

  ADVANCED_SET(${PROJECT_NAME}_VERBOSE_CONFIGURE OFF
    CACHE BOOL
    "Make the ${PROJECT_NAME} configure process verbose."
    )
  
  ADVANCED_SET(${PROJECT_NAME}_ENABLE_EXPLICIT_INSTANTIATION OFF
    CACHE BOOL
    "Enable explicit template instanitation in all packages that support it"
    )
  
  ADVANCED_OPTION(BUILD_SHARED_LIBS "Build shared libraries." OFF)
  
  ADVANCED_SET(${PROJECT_NAME}_INSTALL_INCLUDE_DIR "include"
    CACHE PATH
    "Location where the headers will be installed.  If given as a relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'include'"
    )
  
  ADVANCED_SET(${PROJECT_NAME}_INSTALL_LIB_DIR "lib"
    CACHE PATH
    "Location where the libraries will be installed.  If given as a relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'lib'"
    )
  
  ADVANCED_SET(${PROJECT_NAME}_INSTALL_RUNTIME_DIR "bin"
    CACHE PATH
    "Location where the runtime DLLs will be installed.  If given as a relative path, it will be relative to ${CMAKE_INSTALL_PREFIX}.  If given as an absolute path, it will used as such.  Default is 'bin'"
    )
  
  IF(WIN32 AND NOT CYGWIN)
    SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT OFF)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT ON)
  ENDIF()
  
  ADVANCED_SET(${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES
    ${${PROJECT_NAME}_ENABLE_EXPORT_MAKEFILES_DEFAULT}
    CACHE BOOL
    "Determines if export makefiles will be create and installed."
    )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE OFF CACHE BOOL
    "Allow secondary stable packages and code to be implicitly enabled." )
  
  ADVANCED_SET(${PROJECT_NAME}_TEST_CATEGORIES NIGHTLY CACHE STRING
    "List of categories of tests to enable: '${${PROJECT_NAME}_VALID_CATEGORIES_STR}' (default NIGHLY)."
    )
  PACKAGE_ARCH_ASSERT_VALID_CATEGORIES(${${PROJECT_NAME}_TEST_CATEGORIES})
  
  ADVANCED_SET(${PROJECT_NAME}_REL_CPU_SPEED 1.0 CACHE STRING
    "Relative CPU speed of the computer used to scale performance tests (default 1.0)."
    )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE_DEFAULT}
    CACHE BOOL
    "Determines if a variety of development mode checks are turned on by default or not." )

  ADVANCED_SET( ${PROJECT_NAME}_ASSERT_MISSING_PACKAGES
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL
    "Determines if asserts are performed on missing packages or not." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_STRONG_C_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C code for supported compilers." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_STRONG_CXX_COMPILE_WARNINGS
    ${${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE}
    CACHE BOOL "Enable strong compiler warnings for C++ code for supported compilers." )

  MULTILINE_SET( ENABLE_SHADOW_WARNINGS_DOC
    "Turn ON or OFF shadowing warnings for all packages where strong warnings have"
    " not been explicitly disabled.  Setting the empty '' let's each package decide." )
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS ""
    "${ENABLE_SHADOW_WARNINGS_DOC}" )
  MARK_AS_ADVANCED(${PROJECT_NAME}_ENABLE_SHADOW_WARNINGS)

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_COVERAGE_TESTING OFF
    CACHE BOOL "Enable support for coverage testing by setting needed compiler/linker options." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_CHECKED_STL OFF
    CACHE BOOL "Turn on checked STL checking (e.g. -D_GLIBCXX_DEBUG) or not." )

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_DEBUG_SYMBOLS OFF
    CACHE BOOL "Turn on debugging symbols (e.g. -g) or not if not a full debug build." )

  IF (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "-Werror")
  ELSE()
    SET(${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT "")
  ENDIF()

  ADVANCED_SET( ${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS
    "${${PROJECT_NAME}_WARNINGS_AS_ERRORS_FLAGS_DEFAULT}"
    CACHE STRING "Flags for treating warnings as errors (for all compilers, -Werror by default for GNU).  To turn off warnings as errors set to ''")

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_CIRCULAR_REF_DETECTION_FAILURE OFF CACHE BOOL
    "If test output complaining about circular references is found, then the test will fail." )

  IF (WIN32 AND NOT CYGWIN)
    SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES_DEFAULT FALSE)
  ELSE()
    SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES_DEFAULT TRUE)
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES
    "${${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES_DEFAULT}"
    CACHE BOOL
    "Output any XML dependency files or not." )

  # 2009/01/19: rabartl: Above: This file outputs just fine on MS Windows
  # using MS Visual Studio but it causes the entire file to
  # diff.  There must be something wrong with a newlines or something
  # that is causing this.  If people are going to be doing real
  # development work on MS Windows with MS Visual Studio, then we need
  # to fix this so that the dependency files will get created and
  # checked in correctly.  I will look into this later.

  ADVANCED_SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/python/data/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}"
    CACHE STRING
    "Output XML file containing ${PROJECT_NAME} dependenices used by tools (if not empty)." )
  
  IF(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE)
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/python/data/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
  ELSE()
    SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT "")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
    "${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "Output XML file used by CDash in ${PROJECT_NAME}-independent format (if not empty)." )
  
  IF(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE AND PYTHON_EXECUTABLE)
    SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT
      "${CMAKE_CURRENT_SOURCE_DIR}/cmake/python/data/${${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME}" )
  ELSE()
    SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT "")
  ENDIF()
  ADVANCED_SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
    "${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE_DEFAULT}"
    CACHE STRING
    "HTML ${PROJECT_NAME} dependenices file that will be written to (if not empty)." )

  ADVANCED_SET(${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR
    "" CACHE PATH
    "Output the full XML dependency files in the given directory." )

  ADVANCED_SET(${PROJECT_NAME}_ENABLE_KNOWN_EXTERNAL_REPOS_TYPE
    ""
    CACHE STRING
    "Type of testing to pull in extra respositories (Continuous, or Nightly)" )

  ADVANCED_SET(${PROJECT_NAME}_EXTRAREPOS_FILE
    "${${PROJECT_NAME}_DEPS_HOME_DIR}/cmake/${${PROJECT_NAME}_EXTRA_EXTERNAL_REPOS_FILE_NAME}"
    CACHE FILENAME
    "File contining the list of extra repositories contining add-on packages to process")

  ADVANCED_SET(${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES
    FALSE CACHE BOOL
   "Set if to ignore missing extra repositories (or fail hard)" )

  MESSAGE("")
  MESSAGE("Reading the list of extra repositories from ${${PROJECT_NAME}_EXTRAREPOS_FILE}")
  MESSAGE("")

  INCLUDE(${${PROJECT_NAME}_EXTRAREPOS_FILE})
  PACKAGE_ARCH_PROCESS_EXTRAREPOS_LISTS() # Sets ${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT

  ADVANCED_SET(${PROJECT_NAME}_EXTRA_REPOSITORIES
    "${${PROJECT_NAME}_EXTRA_REPOSITORIES_DEFAULT}"
    CACHE STRING
    "List of external repositories that contain extra ${PROJECT_NAME} packages."
    )

  ADVANCED_SET(${PROJECT_NAME}_INSTALLATION_DIR
    ""
    CACHE STRING
    "Location of an installed version of ${PROJECT_NAME} that will be built against during installation testing"
    )

  IF("${${PROJECT_NAME}_INSTALLATION_DIR}" STREQUAL "")
    SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT OFF)
  ELSE()
    SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT ON)
  ENDIF()
  
  ADVANCED_SET(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING
    ${${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING_DEFAULT}
    CACHE STRING
    "Enable testing against an installed version of ${PROJECT_NAME}."
    )
  
  MARK_AS_ADVANCED(BUILD_TESTING)
  MARK_AS_ADVANCED(CMAKE_BACKWARDS_COMPATIBILITY)
  MARK_AS_ADVANCED(DART_TESTING_TIMEOUT)
  MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH)
  MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH)
  MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES)
  MARK_AS_ADVANCED(CMAKE_OSX_SYSROOT)

ENDMACRO()


#
# Macro that processes the list of TPLs
#

MACRO(PACKAGE_ARCH_PROCESS_TPLS_LISTS)

  SET(TPL_NAME_OFFSET 0)
  SET(TPL_FINDMOD_OFFSET 1)
  SET(TPL_CLASSIFICATION_OFFSET 2)
  SET(TPL_NUM_COLUMNS 3)

  CMAKE_POLICY(PUSH)
  CMAKE_POLICY(SET CMP0007 NEW) # Count empty array elements

  LIST(LENGTH ${PROJECT_NAME}_TPLS_FINDMODS_CLASSIFICATIONS
    ${PROJECT_NAME}_CURR_NUM_TPLS_FULL)
  MATH(EXPR ${PROJECT_NAME}_CURR_NUM_TPLS
    "${${PROJECT_NAME}_CURR_NUM_TPLS_FULL}/${TPL_NUM_COLUMNS}")

  IF (${PROJECT_NAME}_CURR_NUM_TPLS GREATER 0)

    MATH(EXPR ${PROJECT_NAME}_LAST_TPL_IDX
      "${${PROJECT_NAME}_CURR_NUM_TPLS}-1")
    
    IF (NOT APPEND_TO_TPLS_LIST)
      SET(${PROJECT_NAME}_TPLS)
    ENDIF()
  
    FOREACH(TPL_IDX RANGE ${${PROJECT_NAME}_LAST_TPL_IDX})
  
      # Get fields for this TPL
  
      MATH(EXPR TPL_NAME_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_NAME_OFFSET}")
      LIST(GET ${PROJECT_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_NAME_IDX}
        TPL_NAME)
  
      MATH(EXPR TPL_FINDMOD_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_FINDMOD_OFFSET}")
      LIST(GET ${PROJECT_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_FINDMOD_IDX}
        TPL_FINDMOD)
  
      MATH(EXPR TPL_CLASSIFICATION_IDX
        "${TPL_IDX}*${TPL_NUM_COLUMNS}+${TPL_CLASSIFICATION_OFFSET}")
      LIST(GET ${PROJECT_NAME}_TPLS_FINDMODS_CLASSIFICATIONS ${TPL_CLASSIFICATION_IDX}
        TPL_CLASSIFICATION)
  
      # Update TPLS list
  
      LIST(APPEND ${PROJECT_NAME}_TPLS ${TPL_NAME})
  
      # Set ${TPL_NAME}_CLASSIFICATION
  
      IF (TPL_CLASSIFICATION STREQUAL PS
        OR TPL_CLASSIFICATION STREQUAL SS
        OR TPL_CLASSIFICATION STREQUAL TS
        OR TPL_CLASSIFICATION STREQUAL EX
        )
      ELSE()
        MESSAGE(FATAL_ERROR "Error the TPL classification '${TPL_CLASSIFICATION}'"
          " for the TPL ${TPL_NAME} is not a valid classification." )
      ENDIF()
  
      IF (NOT ${TPL_NAME}_CLASSIFICATION) # Allow for testing override
        SET(${TPL_NAME}_CLASSIFICATION ${TPL_CLASSIFICATION})
      ENDIF()
  
      # Set ${TPL_NAME}_FINDMOD
  
      IF (TPL_FINDMOD)
        SET(${TPL_NAME}_FINDMOD ${TPL_FINDMOD})
      ELSE()
        SET(${TPL_NAME}_FINDMOD FindTPL${TPL_NAME})
      ENDIF()
  
      # Set the enable cache variable for ${TPL_NAME}
  
      MULTILINE_SET(DOCSTR
        "Enable support for the TPL ${TPL_NAME} in all supported ${PROJECT_NAME} packages."
        "  This can be set to 'ON', 'OFF', or left empty ''."
        )
      SET_CACHE_ON_OFF_EMPTY( TPL_ENABLE_${TPL_NAME} "" ${DOCSTR} )
  
      # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of
      # ${PROJECT_NAME}_ in order to make it clear that external TPLs are
      # different from packages so users don't get confused and
      # think that the project actually includes some TPL when it does not!
  
    ENDFOREACH()

  ENDIF()
  
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_TPLS)
  ENDIF()

  # Get the final length

  LIST(LENGTH ${PROJECT_NAME}_TPLS ${PROJECT_NAME}_NUM_TPLS)
  PRINT_VAR(${PROJECT_NAME}_NUM_TPLS)
  
  # Create a reverse list for later use
  
  SET(${PROJECT_NAME}_REVERSE_TPLS ${${PROJECT_NAME}_TPLS})
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_TPLS)

  CMAKE_POLICY(POP)

ENDMACRO()


#
# Private helper stuff
#


FUNCTION(PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING PACKAGE_NAME LIST_TYPE
  XML_VAR
  )

  SET(LOC_XML "${${XML_VAR}}")

  SET(DEPS_VAR ${PACKAGE_NAME}_${LIST_TYPE})
  ASSERT_DEFINED(DEPS_VAR)
  SET(DEPS ${${DEPS_VAR}})

  #PRINT_VAR(PACKAGE_NAME)
  #PRINT_VAR(DEPS)

  IF (NOT DEPS)

    APPEND_STRING_VAR(LOC_XML
      "    <${LIST_TYPE}/>\n" )
    
  ELSE()

    SET(VALUE_STR "")

    FOREACH(DEP ${DEPS})

      IF(VALUE_STR)
        SET(VALUE_STR "${VALUE_STR},")
      ENDIF()

      SET(VALUE_STR "${VALUE_STR}${DEP}")

    ENDFOREACH()

    APPEND_STRING_VAR(LOC_XML
      "    <${LIST_TYPE} value=\"${VALUE_STR}\"/>\n" )

  ENDIF()

  IF (LOC_XML)
    SET(${XML_VAR} "${LOC_XML}" PARENT_SCOPE)
  ENDIF()

ENDFUNCTION()


#
# Function that writes the dependency information for ${PROJECT_NAME} into
# an XML file for other tools to use.
#

FUNCTION(PACKAGE_ARCH_DUMP_DEPS_XML_FILE)

  SET(DEPS_XML "")

  APPEND_STRING_VAR(DEPS_XML
    "<PackageDependencies>\n")

  FOREACH(PACKAGE_IDX RANGE ${${PROJECT_NAME}_LAST_PACKAGE_IDX})

    LIST(GET ${PROJECT_NAME}_PACKAGES ${PACKAGE_IDX} PACKAGE)
    LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    
    APPEND_STRING_VAR(DEPS_XML
      "  <Package name=\"${PACKAGE}\" dir=\"${PACKAGE_DIR}\">\n")

    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} LIB_REQUIRED_DEP_PACKAGES DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} LIB_OPTIONAL_DEP_PACKAGES DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} TEST_REQUIRED_DEP_PACKAGES DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} TEST_OPTIONAL_DEP_PACKAGES DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} LIB_REQUIRED_DEP_TPLS DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} LIB_OPTIONAL_DEP_TPLS DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} TEST_REQUIRED_DEP_TPLS DEPS_XML)
    PACKAGE_ARCH_WRITE_DEPS_TO_XML_STRING(${PACKAGE} TEST_OPTIONAL_DEP_TPLS DEPS_XML)

    APPEND_STRING_VAR(DEPS_XML
      "    <EmailAddresses>\n"
      "      <Checkin address=\"${${PACKAGE}_CHECKINS_EMAIL_LIST}\"/>\n"
      "      <Regression address=\"${${PACKAGE}_REGRESSION_EMAIL_LIST}\"/>\n"
      "    </EmailAddresses>\n"
      )

    APPEND_STRING_VAR(DEPS_XML
      "  </Package>\n" )

  ENDFOREACH()

  APPEND_STRING_VAR(DEPS_XML
    "</PackageDependencies>\n" )

  #PRINT_VAR(DEPS_XML)

  FILE(WRITE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ${DEPS_XML} )

ENDFUNCTION()


#
# Macro that ouptuts XML dependency files
#

MACRO(PACKAGE_ARCH_WRITE_XML_DEPENDENCY_FILES)
  
  #PRINT_VAR(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
  IF (${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    IF (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
      SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE})
    ENDIF()
    MESSAGE("" )
    MESSAGE("Dumping the XML dependencies file ${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE} ..." )
    PACKAGE_ARCH_DUMP_DEPS_XML_FILE()
  ENDIF()
  
  #PRINT_VAR(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE)
  IF (${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    IF (NOT IS_ABSOLUTE ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
      SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
        ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE})
    ENDIF()
    MESSAGE("" )
    MESSAGE("Dumping the HTML dependencies webpage file ${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} ..." )
    EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE}
        ${PROJECT_HOME_DIR}/cmake/python/dump-package-dep-table.py
        --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
        --output-html-deps-file=${${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE} )
  ENDIF()
  
  #PRINT_VAR(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE)
  IF (${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE AND ${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE)
    IF (NOT IS_ABSOLUTE ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
      SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE})
    ENDIF()
    MESSAGE("" )
    MESSAGE("Dumping the CDash XML dependencies file ${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE} ..." )
    EXECUTE_PROCESS(
      COMMAND ${PYTHON_EXECUTABLE}
        ${PROJECT_HOME_DIR}/cmake/python/dump-cdash-deps-xml-file.py
        --input-xml-deps-file=${${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE}
        --output-cdash-deps-xml-file=${${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE} )
  ENDIF()

ENDMACRO()


#
# Read in ${PROJECT_NAME} packages and TPLs, process dependencies, write XML files
#
# The reason that these steps are all jammed into one macro is so that the XML
# dependencies of just the core ${PROJECT_NAME} packages can be processed, have the
# XML files written, and then read in the extra set of packages and process
# the dependencies again.
#

MACRO(PACKAGE_ARCH_READ_PACKAGES_PROCESS_DEPENDENCIES_WRITE_XML)

  #
  # 1) Define the lists of all ${PROJECT_NAME} packages and TPLs
  #
  
  # 1.a) Read the core ${PROJECT_NAME} packages

  SET(${PROJECT_NAME}_PACKAGES_FILE
    "${${PROJECT_NAME}_DEPS_HOME_DIR}/cmake/${${PROJECT_NAME}_PACKAGES_FILE_NAME}")

  MESSAGE("")
  MESSAGE("Reading the list of packages from ${${PROJECT_NAME}_PACKAGES_FILE}")
  MESSAGE("")
  
  INCLUDE(${${PROJECT_NAME}_PACKAGES_FILE})
  
  SET(APPEND_TO_PACKAGES_LIST FALSE)
  PACKAGE_ARCH_PROCESS_PACKAGES_AND_DIRS_LISTS()
  
  # 1.b) Read the core TPLs dependencies

  SET(${PROJECT_NAME}_TPLS_FILE
    "${${PROJECT_NAME}_DEPS_HOME_DIR}/cmake/${${PROJECT_NAME}_TPLS_FILE_NAME}")
  
  MESSAGE("")
  MESSAGE("Reading the list of TPLs from ${${PROJECT_NAME}_TPLS_FILE}")
  MESSAGE("")
  
  INCLUDE(${${PROJECT_NAME}_TPLS_FILE})
  
  PACKAGE_ARCH_PROCESS_TPLS_LISTS()
    
  #
  # 2) Process the package and TPL dependencies
  #
  
  PACKAGE_ARCH_READ_ALL_PACKAGE_DEPENDENCIES()

  #
  # 3) Write the XML dependency files for the core ${PROJECT_NAME} packages
  #
  
  IF (${PROJECT_NAME}_OUTPUT_DEPENDENCY_FILES)
    PACKAGE_ARCH_WRITE_XML_DEPENDENCY_FILES()
  ENDIF()

  #
  # 4) Read in the list of externally defined packages and TPLs in external
  # repositories
  #

  # Allow list to be seprated by ',' instead of just by ';'.  This is needed
  # by the unit test driver code
  SPLIT("${${PROJECT_NAME}_EXTRA_REPOSITORIES}"  "," ${PROJECT_NAME}_EXTRA_REPOSITORIES)

  FOREACH(EXTRA_REPO ${${PROJECT_NAME}_EXTRA_REPOSITORIES})

    #PRINT_VAR(EXTRA_REPO)

    # Read in the add-on packages from the extra repo

    SET(EXTRAREPO_PACKAGES_FILE
      "${${PROJECT_NAME}_DEPS_HOME_DIR}/${EXTRA_REPO}/${${PROJECT_NAME}_EXTRA_PACKAGES_FILE_NAME}")

    MESSAGE("")
    MESSAGE("Reading a list of extra packages from ${EXTRAREPO_PACKAGES_FILE} ... ")
    MESSAGE("")

    IF (NOT EXISTS "${EXTRAREPO_PACKAGES_FILE}" AND ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
      MESSAGE(
        "\n***"
        "\n** WARNING!  Ignoring missing extra repo '${EXTRAREPO_PACKAGES_FILE}' on request!"
        "\n***\n")
    ELSE()
      INCLUDE("${EXTRAREPO_PACKAGES_FILE}")  # Writes the variable ???
      SET(APPEND_TO_PACKAGES_LIST TRUE)
      PACKAGE_ARCH_PROCESS_PACKAGES_AND_DIRS_LISTS()  # Reads the variable ???
    ENDIF()

    # Read in the add-on TPLs from the extra repo

    SET(EXTRAREPO_TPLS_FILE
      "${${PROJECT_NAME}_DEPS_HOME_DIR}/${EXTRA_REPO}/${${PROJECT_NAME}_EXTRA_TPLS_FILE_NAME}")

    MESSAGE("")
    MESSAGE("Reading a list of extra TPLs from ${EXTRAREPO_TPLS_FILE} ... ")
    MESSAGE("")

    IF (NOT EXISTS "${EXTRAREPO_TPLS_FILE}" AND ${PROJECT_NAME}_IGNORE_MISSING_EXTRA_REPOSITORIES)
      MESSAGE(
        "\n***"
        "\n** WARNING!  Ignoring missing extra repo '${EXTRAREPO_TPLS_FILE}' on request!"
        "\n***\n")
    ELSE()
      INCLUDE("${EXTRAREPO_TPLS_FILE}")  # Writes the varaible ???
      SET(APPEND_TO_TPLS_LIST TRUE)
      PACKAGE_ARCH_PROCESS_TPLS_LISTS()  # Reads the variable ???
    ENDIF()

  ENDFOREACH()

  #
  # 5) Read in the package dependencies again to now pick up all of the
  # defined packages (not just the core packages)
  #

  IF (${PROJECT_NAME}_EXTRA_REPOSITORIES)
    PACKAGE_ARCH_READ_ALL_PACKAGE_DEPENDENCIES()
  ENDIF()

  #
  # 6) Write out the XML dependency files again but this time for the full
  # list in the build directory!
  #

  IF (${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR)
    #MESSAGE("Writing dependency XML files in ${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR} ...")
    SET(${PROJECT_NAME}_DEPS_XML_OUTPUT_FILE
      "${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_XML_FILE_NAME}" )
    IF(PYTHON_EXECUTABLE)
      SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE
        "${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}/${${PROJECT_NAME}_CDASH_SUBPROJECT_DEPS_XML_FILE_NAME}" )
      SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE
        "${${PROJECT_NAME}_OUTPUT_FULL_DEPENDENCY_FILES_IN_DIR}/${${PROJECT_NAME}_PACKAGE_DEPS_TABLE_HTML_FILE_NAME}" )
    ELSE()
      SET(${PROJECT_NAME}_CDASH_DEPS_XML_OUTPUT_FILE "")
      SET(${PROJECT_NAME}_DEPS_HTML_OUTPUT_FILE "")
    ENDIF()
    PACKAGE_ARCH_WRITE_XML_DEPENDENCY_FILES()
  ENDIF()

ENDMACRO()


#
# Adjust package enable logic and print out before and after state
#
# On output sets:
#
#    ${PROJECT_NAME}_NUM_ENABLED_PACKAGES: Number of enabled packages (local variable)
#    ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}: Enable status of PACKAGE_NAME (local variable)
#    ToDo: Fill in others as well!
#

MACRO(PACKAGE_ARCH_ADJUST_AND_PRINT_PACKAGE_DEPENDENCIES)

  PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
    "\nExplicitly enabled packages on input (by user)" ON FALSE)
  PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
    "\nExplicitly disabled packages on input (by user or by default)" OFF FALSE)
  PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
    "\nExplicitly enabled TPLs on input (by user)" ON FALSE)
  PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
    "\nExplicitly disabled TPLs on input (by user or by default)" OFF FALSE)

  SET(DO_PROCESS_MPI_ENABLES ON) # Should not be needed but CMake is not working!
  PACKAGE_ARCH_ADJUST_PACKAGE_ENABLES(TRUE)

  PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
    "\nFinal set of enabled packages" ON FALSE)
  PACKAGE_ARCH_PRINT_ENABLED_PACKAGE_LIST(
    "\nFinal set of non-enabled packages" OFF TRUE)
  PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
    "\nFinal set of enabled TPLs" ON FALSE)
  PACKAGE_ARCH_PRINT_ENABLED_TPL_LIST(
    "\nFinal set of non-enabled TPLs" OFF TRUE)

  PACKAGE_ARCH_GET_ENABLED_LIST( ${PROJECT_NAME}_PACKAGES ${PROJECT_NAME} ON FALSE
    ${PROJECT_NAME}_ENABLED_PACKAGES ${PROJECT_NAME}_NUM_ENABLED_PACKAGES )

ENDMACRO()


#
# Macro that gathers information from enabled TPLs
#

MACRO(PACKAGE_ARCH_PROCESS_ENABLED_TPLS)
  FOREACH(TPL_NAME ${${PROJECT_NAME}_TPLS})
    IF (TPL_ENABLE_${TPL_NAME})
      MESSAGE(STATUS "Processing enabled TPL: ${TPL_NAME}")
      INCLUDE(TPLs/${${TPL_NAME}_FINDMOD})
      ASSERT_DEFINED(TPL_${TPL_NAME}_INCLUDE_DIRS)
      ASSERT_DEFINED(TPL_${TPL_NAME}_LIBRARIES)
      ASSERT_DEFINED(TPL_${TPL_NAME}_LIBRARY_DIRS)
    ENDIF()
  ENDFOREACH()
ENDMACRO()


#
# Macros for setting up the standard environment
#


#
# Pre-setup of the environment
#

MACRO(PACKAGE_ARCH_PRE_SETUP_ENV)

  # Set to release build by default
  
  IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting CMAKE_BUILD_TYPE=RELEASE since it was not set ...")
    SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Type of build to perform (i.e. DEBUG, RELEASE, NONE)" )
  ELSE()
    STRING(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UP)
    LIST(FIND CMAKE_BUILD_TYPES_LIST ${CMAKE_BUILD_TYPE_UP} BUILD_TYPE_IDX)
    IF (BUILD_TYPE_IDX EQUAL -1)
      MESSAGE(SEND_ERROR "Error, the given CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
        " is not in the list of valid values \"${CMAKE_BUILD_TYPES_LIST}\"!")
    ENDIF()
  ENDIF()
  PRINT_VAR(CMAKE_BUILD_TYPE)

  # Set up MPI if MPI is being used

  ASSERT_DEFINED(TPL_ENABLE_MPI)
  IF (TPL_ENABLE_MPI)
    PACKAGE_ARCH_SETUP_MPI()
  ENDIF()

  # Enable compilers
  
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_C)
  IF (${PROJECT_NAME}_ENABLE_C)
    ENABLE_LANGUAGE(C)
    INCLUDE(CMakeDetermineCCompiler)
    PRINT_VAR(CMAKE_C_COMPILER_ID)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  ENDIF()
  
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX)
  IF (${PROJECT_NAME}_ENABLE_CXX)
    ENABLE_LANGUAGE(CXX)
    INCLUDE(CMakeDetermineCXXCompiler)
    PRINT_VAR(CMAKE_CXX_COMPILER_ID)
    # See CMake/Modules/CMakeCXXCompilerId.cpp.in in the CMake source
    # directory for a listing of known compiler types.
  ENDIF()
  
  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_Fortran)
  IF (${PROJECT_NAME}_ENABLE_Fortran)
    ENABLE_LANGUAGE(Fortran)
  ENDIF()

  # Set up for strong compiler warnings and warnings as errors
 
  INCLUDE(PackageArchSetupBasicCompileLinkFlags)
  PACKAGE_ARCH_SETUP_BASIC_COMPILE_LINK_FLAGS()

  # Find the host site name used in selecting or deselecting tests by the
  # PACKAGE_ADD_TEST(...) function.
  
  SITE_NAME(${PROJECT_NAME}_HOSTNAME)
  MARK_AS_ADVANCED(${PROJECT_NAME}_HOSTNAME)
  PRINT_VAR(${PROJECT_NAME}_HOSTNAME)

  # Find the host site type name used in selecting or deselecting tests by the
  # PACKAGE_ADD_TEST(...) function.

  PRINT_VAR(CMAKE_HOST_SYSTEM_NAME)

ENDMACRO()


#
# Post-setup of the environment
#

MACRO(PACKAGE_ARCH_POST_SETUP_ENV)

  # Set the hack library to get link options on

  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Creating dummy last_lib for appending the link flags: "
        "${${PROJECT_NAME}_EXTRA_LINK_FLAGS}")
    ENDIF()
    IF (NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c)
      FILE(WRITE ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c
        "typedef int last_lib_dummy_t;\n")
    ENDIF()
    ADD_LIBRARY(last_lib STATIC ${CMAKE_CURRENT_BINARY_DIR}/last_lib_dummy.c)
    TARGET_LINK_LIBRARIES(last_lib ${${PROJECT_NAME}_EXTRA_LINK_FLAGS})
  ENDIF()

ENDMACRO()


#
# Macro that does the final set of package configurations
#

MACRO(PACKAGE_ARCH_CONFIGURE_ENABLED_PACKAGES)

  GLOBAL_NULL_SET(${PROJECT_NAME}_INCLUDE_DIRS)
  GLOBAL_NULL_SET(${PROJECT_NAME}_LIBRARY_DIRS)
  GLOBAL_NULL_SET(${PROJECT_NAME}_LIBRARIES)

  SET(CONFIGURED_A_PACKAGE FALSE)
  SET(ENABLED_PACKAGE_LIBS_TARGETS)
  
  FOREACH(PACKAGE_IDX RANGE ${${PROJECT_NAME}_LAST_PACKAGE_IDX})
    LIST(GET ${PROJECT_NAME}_PACKAGES ${PACKAGE_IDX} PACKAGE)
    LIST(GET ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    IF(${PROJECT_NAME}_ENABLE_${PACKAGE})
      SET(PACKAGE_NAME_GLOBAL ${PACKAGE}) # For consistency checking
      SET(${PACKAGE}_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/packages/${PACKAGE_DIR})
      IF(${PACKAGE}_SPECIFIED_BUILD_DIR)
        IF(IS_ABSOLUTE ${${PACKAGE}_SPECIFIED_BUILD_DIR})
          SET(${PACKAGE}_BINARY_DIR ${${PACKAGE}_SPECIFIED_BUILD_DIR})
        ELSE()
          SET(${PACKAGE}_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/${${PACKAGE}_SPECIFIED_BUILD_DIR})
        ENDIF()
      ELSE()
	SET(${PACKAGE}_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/packages/${PACKAGE_DIR})
      ENDIF()
      IF (NOT EXISTS ${${PACKAGE}_SOURCE_DIR}/CMakeLists.txt)
        MESSAGE(FATAL_ERROR "Error, the file ${PACKAGE_DIR}/CMakeLists.txt does not exist!")
      ENDIF()
      ADD_SUBDIRECTORY(${${PACKAGE}_SOURCE_DIR} ${${PACKAGE}_BINARY_DIR})
      LIST(APPEND ENABLED_PACKAGE_LIBS_TARGETS ${PACKAGE}_libs)
      LIST(APPEND ${PROJECT_NAME}_INCLUDE_DIRS ${${PACKAGE}_INCLUDE_DIRS})
      LIST(APPEND ${PROJECT_NAME}_LIBRARY_DIRS ${${PACKAGE}_LIBRARY_DIRS})
      LIST(APPEND ${PROJECT_NAME}_LIBRARIES ${${PACKAGE}_LIBRARIES})
      SET(CONFIGURED_A_PACKAGE TRUE)
    ENDIF()
  ENDFOREACH()

  ADVANCED_SET( ${PROJECT_NAME}_ALLOW_NO_PACKAGES ON
    CACHE BOOL "Allow configuration to finish even if no packages are enabled")

  IF (NOT CONFIGURED_A_PACKAGE)
    IF (${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      SET(MSG_TYPE WARNING)
    ELSE()
      SET(MSG_TYPE ERROR)
    ENDIF()
    MESSAGE(
      "\n***"
      "\n*** ${MSG_TYPE}:  There were no packages configured so no libraries"
        " or tests/examples will be built!"
      "\n***\n"
      )
    IF (NOT ${PROJECT_NAME}_ALLOW_NO_PACKAGES)
      MESSAGE(SEND_ERROR "Stopping configure!")
    ENDIF()
  ELSE()
    ASSERT_AND_TOUCH_DEFINED(${PROJECT_NAME}_ALLOW_NO_PACKAGES)
  ENDIF()
  
  REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_INCLUDE_DIRS)
  REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_LIBRARY_DIRS)
  REMOVE_GLOBAL_DUPLICATES(${PROJECT_NAME}_LIBRARIES)

  # Add global 'libs;' target
  IF(NOT ENABLED_PACKAGE_LIBS_TARGETS)
    RETURN()
  ENDIF()
  LIST(REVERSE ENABLED_PACKAGE_LIBS_TARGETS)
  # make it so when no packages are enabled 
  # it is not a cmake error
  IF (${PROJECT_NAME}_EXTRA_LINK_FLAGS)
    APPEND_SET(ENABLED_PACKAGE_LIBS_TARGETS last_lib)
  ENDIF()
  #PRINT_VAR(ENABLED_PACKAGE_LIBS_TARGETS)
  ADD_CUSTOM_TARGET(${PROJECT_NAME}_libs)
  ADD_DEPENDENCIES(${PROJECT_NAME}_libs ${ENABLED_PACKAGE_LIBS_TARGETS})
  ADD_CUSTOM_TARGET(libs)
  ADD_DEPENDENCIES(libs ${ENABLED_PACKAGE_LIBS_TARGETS})

ENDMACRO()


#
#  Macro that allows packages to easily make a feature SS for development
#  builds and PS for release builds
#.
#  The OUTPUT_VAR is set to ON or OFF based on the configure state. In
#  development mode it will be set to ON only if SS code is enabled, 
#  otherwise it is set to OFF. In release mode it is always set to ON.
#  This allows some sections of PROJECT_NAME to be considered SS for 
#  development mode reducing testing time, while still having important
#  functionality available to users by default
MACRO(PACKAGE_ARCH_SET_SS_FOR_DEV_MODE OUTPUT_VAR)
  IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
    SET(${OUTPUT_VAR} ${${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE})
  ELSE()
    SET(${OUTPUT_VAR} ON)
  ENDIF()
ENDMACRO()


#
# Macro that drives a experimental 'dashboard' target
#

MACRO(PACKAGE_ARCH_ADD_DASHBOARD_TARGET)

  IF (NOT (WIN32 AND NOT CYGWIN))

    ADVANCED_SET(${PROJECT_NAME}_DASHBOARD_CTEST_ARGS "" CACHE STRING
      "Extra arguments to pass to CTest when calling 'ctest -S' to run the 'dashboard' make target." )

    ADVANCED_SET(CTEST_BUILD_FLAGS "" CACHE STRING
      "Sets CTEST_BUILD_FLAGS on the env before invoking 'ctest -S'." )

    ADVANCED_SET(CTEST_PARALLEL_LEVEL "" CACHE STRING
      "Sets CTEST_PARALLEL_LEVEL on the env before invoking 'ctest -S'." )
  
    # H.1) Enable all packages that are enabled and have tests enabled
  
    SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
    SET(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST)
    FOREACH(PACKAGE ${${PROJECT_NAME}_PACKAGES})
      IF (${PROJECT_NAME}_ENABLE_${PACKAGE} AND ${PACKAGE}_ENABLE_TESTS)
        IF (${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
          SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST
            "${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}\;${PACKAGE}") 
        ELSE()
          SET(${PROJECT_NAME}_ENABLED_PACKAGES_LIST "${PACKAGE}") 
        ENDIF()
        SET(${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST
          ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST} -D${PROJECT_NAME}_ENABLE_${PACKAGE}=ON)
      ENDIF()
    ENDFOREACH()
    #PRINT_VAR(${PROJECT_NAME}_ENABLED_PACKAGES_LIST)
    
    SET(EXPR_CMND_ARGS)
    IF (${PROJECT_NAME}_ENABLE_COVERAGE_TESTING)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DO_COVERAGE_TESTING=TRUE")
    ENDIF()
    IF (CTEST_BUILD_FLAGS)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_BUILD_FLAGS='${CTEST_BUILD_FLAGS}'")
    ENDIF()
    IF (CTEST_PARALLEL_LEVEL)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_PARALLEL_LEVEL='${CTEST_PARALLEL_LEVEL}'")
    ENDIF()
    IF (CTEST_DROP_SITE)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_SITE=${CTEST_DROP_SITE}")
    ENDIF()
    IF (CTEST_DROP_LOCATION)
      APPEND_SET(EXPR_CMND_ARGS "CTEST_DROP_LOCATION=${CTEST_DROP_LOCATION}")
    ENDIF()

    #PRINT_VAR(${PROJECT_NAME}_EXTRA_REPOSITORIES)
    JOIN(${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED "," FALSE
      ${${PROJECT_NAME}_EXTRA_REPOSITORIES})
    APPEND_SET(EXPR_CMND_ARGS
      ${PROJECT_NAME}_EXTRA_REPOSITORIES=${${PROJECT_NAME}_EXTRA_REPOSITORIES_JOINED})
  
    # H.2) Add the custom target to enable all the packages with tests enabled
    
    ADD_CUSTOM_TARGET(dashboard
  
      VERBATIM
    
      # WARNING: The echoed command and the actual commands are duplicated!  You have to reproduce them!
  
      COMMAND echo
      COMMAND echo "***************************************************"
      COMMAND echo "*** Running incremental experimental dashboard ***" 
      COMMAND echo "***************************************************"
      COMMAND echo
      COMMAND echo ${PROJECT_NAME}_ENABLED_PACKAGES_LIST=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
      COMMAND echo
  
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** A) Clean out the list of packages"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_HOME_DIR}
  
      # NOTE: Above, if ${PROJECT_NAME}_ENABLE_ALL_PACKAGES was set in CMakeCache.txt, then setting
      # -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF will turn it off in the cache.  Note that it will
      # never be turned on again which means that the list of packages will be set explicitly below.
    
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** B) Run the dashboard command setting the list of packages"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: env ${EXPR_CMND_ARGS}
        ${PROJECT_NAME}_PACKAGES=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
        ${CMAKE_CTEST_COMMAND} ${${PROJECT_NAME}_DASHBOARD_CTEST_ARGS} -S
          ${PROJECT_HOME_DIR}/cmake/ctest/experimental_build_test.cmake
      COMMAND echo
      COMMAND env ${EXPR_CMND_ARGS}
        ${PROJECT_NAME}_PACKAGES=${${PROJECT_NAME}_ENABLED_PACKAGES_LIST}
        ${CMAKE_CTEST_COMMAND} ${${PROJECT_NAME}_DASHBOARD_CTEST_ARGS} -S
          ${PROJECT_HOME_DIR}/cmake/ctest/experimental_build_test.cmake || echo
  
      # 2009/07/05: rabartl: Above, I added the ending '|| echo' to always make
      # the command pass so that 'make will not stop and avoid this last command
      # to set back the enabled packages.
  
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** C) Clean out the list of packages again to clean the cache file"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} -D${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES:BOOL=TRUE
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON -D${PROJECT_NAME}_ENABLE_ALL_PACKAGES:BOOL=OFF ${PROJECT_HOME_DIR}
    
      COMMAND echo
      COMMAND echo "***"
      COMMAND echo "*** D) Reconfigure with the original package list"
      COMMAND echo "***"
      COMMAND echo
      COMMAND echo Running: ${CMAKE_COMMAND} ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST}
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON ${PROJECT_HOME_DIR}
      COMMAND echo
      COMMAND ${CMAKE_COMMAND} ${${PROJECT_NAME}_ENABLED_PACKAGES_CMAKE_ARG_LIST}
        -D${PROJECT_NAME}_ALLOW_NO_PACKAGES:BOOL=ON ${PROJECT_HOME_DIR}
  
      COMMAND echo
      COMMAND echo "See the results at http://${CTEST_DROP_SITE}${CTEST_DROP_LOCATION}&display=project\#Experimental"
      COMMAND echo
   
      )
  
  ENDIF()

ENDMACRO()


MACRO(PACKAGE_ARCH_EXCLUDE_FILES)
  SET(FILES_TO_EXCLUDE ${ARGN})
  
  #need to add "/<project source dir>/packages/<package name>/" to each file this is to prevent
  #someone from trying to exclude a file like "readme" and having it 
  #inadvertently exclude a file matching that name in another package.
  SET(MODIFIED_FILES_TO_EXCLUDE "")
  FOREACH(FILE ${FILES_TO_EXCLUDE})
    #Ensure that if the full path was specified for the file that we don't add
    #"/<project source dir>/packages/<package name>/" again.
    SET(MATCH_STRING "${${PROJECT_NAME}_SOURCE_DIR}/packages/${PACKAGE_DIR}")
    STRING(REGEX MATCH ${MATCH_STRING} MATCHED ${FILE} )
    IF(NOT MATCHED)
      LIST(APPEND MODIFIED_FILES_TO_EXCLUDE ${${PROJECT_NAME}_SOURCE_DIR}/packages/${PACKAGE_DIR}/${FILE})
    ELSE()
      LIST(APPEND MODIFIED_FILES_TO_EXCLUDE ${FILE})
    ENDIF()
  ENDFOREACH()
 
#Leaving in for debugging purposes
#  MESSAGE("List of files being excluded for package ${PACKAGE_NAME}")
#  FOREACH(NEW_FILE ${MODIFIED_FILES_TO_EXCLUDE})
#    MESSAGE(${NEW_FILE})
#  ENDFOREACH()

  LIST(APPEND CPACK_SOURCE_IGNORE_FILES ${MODIFIED_FILES_TO_EXCLUDE})
  SET(CPACK_SOURCE_IGNORE_FILES ${CPACK_SOURCE_IGNORE_FILES} PARENT_SCOPE)

ENDMACRO()


#
#  macro for helping set up exclude files only for the packages
#  that will not be supporting autotools.
#  Returns a list of the autotools files that shoudl be excluded for
#  the package.
#
#  example: PACKAGE_APPLY_TO_NO_AUTOTOOLS_PACKAGES("configure.ac" list)
#    assuming that the packages epetra and teuchos are not supporting 
#    autotools anymore then the return value would be:
#    "epetra/configure.ac;teuchos/configure.ac"
#
#

MACRO(PACKAGE_ARCH_EXCLUDE_AUTOTOOLS_FILES) # PACKAGE_NAME LIST_RETURN)
  SET(AUTOTOOLS_FILES 
    configure.ac
    configure
    Makefile.am
    Makefile.in
    .*.m4
    bootstrap
    config/
    )

  FOREACH(FILE ${AUTOTOOLS_FILES})
    PACKAGE_ARCH_EXCLUDE_FILES(${FILE} \(.*/\)*${FILE}) 
  ENDFOREACH()
ENDMACRO()
