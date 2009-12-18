
INCLUDE(PackageArchProcessPackagesAndDirsLists)
INCLUDE(PackageArchAdjustPackageEnables)
INCLUDE(PackageArchSetupMPI)

INCLUDE(AddOptionAndDefine)
INCLUDE(AdvancedOption)
INCLUDE(AdvancedSet)
INCLUDE(AppendStringVar)
INCLUDE(CMakeBuildTypesList)
INCLUDE(FindListElement)
INCLUDE(GlobalNullSet)
INCLUDE(PrintNonemptyVar)
INCLUDE(PrintVar)
INCLUDE(RemoveGlobalDuplicates)


#
# Define all of the standard global package architecture options.
#

MACRO(PACKAGE_ARCH_DEFINE_GLOBAL_OPTIONS)

  SET( ${PROJECT_NAME}_ENABLE_ALL_PACKAGES OFF CACHE BOOL
    "Enable all packages (Primary Stable and perhaps Secondary Stable packages)." )
  
  SET(${PROJECT_NAME}_ENABLE_ALL_OPTIONAL_PACKAGES OFF CACHE BOOL
    "Recursively enable all optional packages for set of enabled packages." )
  
  SET(${PROJECT_NAME}_ENABLE_ALL_FORWARD_DEP_PACKAGES OFF CACHE BOOL
    "Recursively enable all packages that have required or optional dependencies for set of enabled packages." )
  
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_TESTS ""
    "Enable tests in all packages  (set to ON, OFF, or leave empty)." )
  
  SET_CACHE_ON_OFF_EMPTY(${PROJECT_NAME}_ENABLE_EXAMPLES ""
    "Enable examples in all packages  (set to ON, OFF, or leave empty).  If left empty, then this will be set to ON if ${PROJECT_NAME}_ENABLE_TESTS=ON" )
  
  IF (${PROJECT_NAME}_ENABLE_TESTS AND ${PROJECT_NAME}_ENABLE_EXAMPLES STREQUAL "")
    MESSAGE(STATUS "Setting ${PROJECT_NAME}_ENABLE_EXAMPLES=ON because ${PROJECT_NAME}_ENABLE_TESTS=ON")
    SET(${PROJECT_NAME}_ENABLE_EXAMPLES ON)
  ENDIF()

  SET( ${PROJECT_NAME}_UNENABLE_ENABLED_PACKAGES OFF CACHE BOOL
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

  ADVANCED_SET( ${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE
    ON  #NOTE: Change this to 'OFF' in a release branch!
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

  LIST(LENGTH ${PROJECT_NAME}_TPLS_AND_CLASSIFICATIONS ${PROJECT_NAME}_NUM_TPLS_2)
  MATH(EXPR ${PROJECT_NAME}_NUM_TPLS "${${PROJECT_NAME}_NUM_TPLS_2}/2")
  PRINT_VAR(${PROJECT_NAME}_NUM_TPLS)
  MATH(EXPR ${PROJECT_NAME}_LAST_TPL_IDX "${${PROJECT_NAME}_NUM_TPLS}-1")
  
  SET(${PROJECT_NAME}_TPLS)

  FOREACH(TPL_IDX RANGE ${${PROJECT_NAME}_LAST_TPL_IDX})

    MATH(EXPR TPL_NAME_IDX "${TPL_IDX}*2")
    MATH(EXPR TPL_ENABLED_IDX "${TPL_IDX}*2+1")

    LIST(GET ${PROJECT_NAME}_TPLS_AND_CLASSIFICATIONS ${TPL_NAME_IDX} TPL)
    LIST(GET ${PROJECT_NAME}_TPLS_AND_CLASSIFICATIONS ${TPL_ENABLED_IDX} TPL_CLASSIFICATION)

    LIST(APPEND ${PROJECT_NAME}_TPLS ${TPL})

    IF (TPL_CLASSIFICATION STREQUAL PS
      OR TPL_CLASSIFICATION STREQUAL SS
      OR TPL_CLASSIFICATION STREQUAL TS
      OR TPL_CLASSIFICATION STREQUAL EX
      )
    ELSE()
      MESSAGE(FATAL_ERROR "Error the TPL classification '${TPL_CLASSIFICATION}'"
        " for the TPL ${TPL} is not a valid classification." )
    ENDIF()

    IF (NOT ${TPL}_CLASSIFICATION) # Allow for testing override
      SET(${TPL}_CLASSIFICATION ${TPL_CLASSIFICATION})
    ENDIF()

    MULTILINE_SET(DOCSTR
      "Enable support for the TPL ${TPL} in all supported ${PROJECT_NAME} packages."
      "  This can be set to 'ON', 'OFF', or left empty ''."
      )
    SET( TPL_ENABLE_${TPL} "" CACHE STRING ${DOCSTR})

    # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of
    # ${PROJECT_NAME}_ in order to make it clear that external TPLs are
    # different from packages so users don't get confused and
    # think that the project actually includes some TPL when it does not!

  ENDFOREACH()
  
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_TPLS)
  ENDIF()
  
  # Create a reverse list for later use
  
  SET(${PROJECT_NAME}_REVERSE_TPLS ${${PROJECT_NAME}_TPLS})
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_TPLS)

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

    SET(ADDRESS_URL_BASE "software.sandia.gov")

    STRING(TOLOWER "${PACKAGE}" LPACKAGE)
    APPEND_STRING_VAR(DEPS_XML
      "    <EmailAddresses>\n"
      "      <Checkin address=\"${LPACKAGE}-checkins@${ADDRESS_URL_BASE}\"/>\n"
      "      <Regression address=\"${LPACKAGE}-regression@${ADDRESS_URL_BASE}\"/>\n"
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
# Set up the standard environment
#


MACRO(PACKAGE_ARCH_PRE_SETUP_ENV)

  # Set to release build by default
  
  IF (NOT CMAKE_BUILD_TYPE)
    MESSAGE(STATUS "Setting CMAKE_BUILD_TYPE=RELEASE since it was not set ...")
    SET(CMAKE_BUILD_TYPE RELEASE)
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
# Macro that gathers information from enabled TPLs
#

MACRO(PACKAGE_ARCH_PROCESS_ENABLED_TPLS)
  FOREACH(TPL ${${PROJECT_NAME}_TPLS})
    IF (TPL_ENABLE_${TPL})
      MESSAGE(STATUS "Processing enabled TPL: ${TPL}")
      INCLUDE(TPLs/FindTPL${TPL})
      ASSERT_DEFINED(TPL_${TPL}_INCLUDE_DIRS)
      ASSERT_DEFINED(TPL_${TPL}_LIBRARIES)
      ASSERT_DEFINED(TPL_${TPL}_LIBRARY_DIRS)
    ENDIF()
  ENDFOREACH()
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
      IF (NOT EXISTS ${PROJECT_SOURCE_DIR}/packages/${PACKAGE_DIR}/CMakeLists.txt)
        MESSAGE(FATAL_ERROR "Error, the file ${PACKAGE_DIR}/CMakeLists.txt does not exist!")
      ENDIF()
      ADD_SUBDIRECTORY(packages/${PACKAGE_DIR})
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
