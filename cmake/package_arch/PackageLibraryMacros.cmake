INCLUDE(PackageCreateClientTemplateHeaders)
INCLUDE(ParseVariableArguments)
INCLUDE(GlobalSet)
INCLUDE(AppendSet)
INCLUDE(AppendGlob)
INCLUDE(AppendGlobalSet)
INCLUDE(PrependGlobalSet)
INCLUDE(RemoveGlobalDuplicates)
INCLUDE(PackageGeneralMacros)
INCLUDE(SetAndIncDirs)


#
# Macro that configures the package's main config.h file
#

FUNCTION(PACKAGE_CONFIGURE_FILE PACKAGE_NAME_CONFIG_FILE)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  ENDIF()

  STRING(TOUPPER "${PACKAGE_NAME}" PACKAGE_NAME_UC)

  IF (${PACKAGE_NAME}_SHOW_DEPRECATED_WARNINGS)
    MULTILINE_SET(${PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS
      "#ifndef ${PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))\n"
      "#    define ${PACKAGE_NAME_UC}_DEPRECATED  __attribute__((__deprecated__))\n"
      "#  else\n"
      "#    define ${PACKAGE_NAME_UC}_DEPRECATED\n"
      "#  endif\n"
      "#endif\n"
      )
  ELSE()
    SET(${PACKAGE_NAME_UC}_DEPRECATED_DECLARATIONS "#define ${PACKAGE_NAME_UC}_DEPRECATED")
  ENDIF()

  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

ENDFUNCTION()


#
# Macro used to add a package library
#

FUNCTION(PACKAGE_ADD_LIBRARY LIBRARY_NAME)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_LIBRARY: ${LIBRARY_NAME}")
    IF(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
      MESSAGE("\n${PACKAGE_NAME}_LIBRARIES In installation testing mode, libraries will be found instead of created.")
    ENDIF()
  ENDIF()

  PARSE_ARGUMENTS(
    PARSE #prefix
    "HEADERS;NOINSTALLHEADERS;SOURCES;DEPLIBS;IMPORTEDLIBS;DEFINES" # Lists
    "TESTONLY;NO_INSTALL_LIB_OR_HEADERS;CUDALIBRARY" #Options
    ${ARGN} # Remaining arguments passed in
    )

  #if we are doing installation testing we want to skip adding libraries, unless
  #they are test only libraries which are not installed. 
  IF(NOT ${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING OR PARSE_TESTONLY)

    # Add the link directory for this library.

    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS ${CMAKE_CURRENT_BINARY_DIR})

    # NOTE: Above , this link path not really used here for anything.
    # Instead it is just added to the other set link library directories
    # that are already set.  These link directories are then extracted
    # and stored into stored in ${PACKAGE_NAME}_LIBRARY_DIRS.

    # Add whatever include directories have been defined so far

    INCLUDE_DIRECTORIES(AFTER ${${PACKAGE_NAME}_INCLUDE_DIRS})

    # Add whatever link directories have been added so far

    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS ${${PACKAGE_NAME}_LIBRARY_DIRS})
    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS ${${PACKAGE_NAME}_TEST_LIBRARY_DIRS})

    # Add dependent libraries passed directly in

    IF (PARSE_DEPLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "DEPLIBS = ${PARSE_DEPLIBS}")
    ENDIF()
    IF (PARSE_IMPORTEDLIBS AND ${PROJECT_NAME}_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "IMPORTEDLIBS = ${PARSE_IMPORTEDLIBS}")
    ENDIF()

    SET(LINK_LIBS)

    IF (PARSE_DEPLIBS)
      APPEND_SET(LINK_LIBS ${PARSE_DEPLIBS})
    ENDIF()
    IF (PARSE_IMPORTEDLIBS)
      APPEND_SET(LINK_LIBS ${PARSE_IMPORTEDLIBS})
    ENDIF()

    # We only want to link to the dependent package and TPL libraries when we need
    # to.  We only need to link to these dependent libraries when this is the first
    # library being created for this package or if this library does not depend
    # on other libraries created for this package.  Otherwise, we don't need to
    # add the include directories or link libraries because a dependent lib
    # specified in PARSE_DEP_LIBS already has everything that we need.

    # We also Need to make special considerations for test libraries since
    # things need to be handled a little bit differently (but not much).  In the
    # case of test libaries, we need to also pull the test-only dependencies.
    # In this case, we will always assume that we will add in the test
    # libraries.

    SET(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

    IF (PARSE_DEPLIBS AND NOT PARSE_TESTONLY)
      FOREACH(DEPLIB ${PARSE_DEPLIBS})
        LIST(FIND ${PACKAGE_NAME}_LIBRARIES ${DEPLIB} DEPLIB_IDX)
        IF (NOT DEPLIB_IDX EQUAL -1)
          # The library being created here is dependent on another of this package's
          # libraries so there is no need to add in this package's dependent package
          # and TPL libraries.
          SET(ADD_DEP_PACKAGE_AND_TPL_LIBS FALSE)
        ENDIF()
      ENDFOREACH()
    ELSE()
      # If there are no dependent libs passed in, then this library can not possiblly
      # depend on the package's other libraries so we must link to the dependent libraries
      # in dependent libraries and TPLs.
    ENDIF()

    IF (ADD_DEP_PACKAGE_AND_TPL_LIBS)

      IF (NOT PARSE_TESTONLY)
        SET(TEST_OR_LIB_ARG LIB)
      ELSE()
        SET(TEST_OR_LIB_ARG TEST)
      ENDIF()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "\nPulling in header and libraries dependencies for ${TEST_OR_LIB_ARG} ...")
      ENDIF()

      # Add the dependent package libraries (if we have not done so yet for this package)
      PACKAGE_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} ${TEST_OR_LIB_ARG} PACKAGES ALL_DEP_PACKAGES)
      PACKAGE_SORT_AND_APPEND_PATHS_LIBS("${${PROJECT_NAME}_REVERSE_PACKAGES}"
        "${ALL_DEP_PACKAGES}" "" LINK_LIBS)

      # Add the TPL libraries (if we have not done so yet for this package)
      PACKAGE_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} ${TEST_OR_LIB_ARG} TPLS ALL_TPLS)
      PACKAGE_SORT_AND_APPEND_PATHS_LIBS("${${PROJECT_NAME}_REVERSE_TPLS}"
       "${ALL_TPLS}" TPL_ LINK_LIBS)

    ENDIF()

    IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
      PRINT_VAR(LINK_LIBS)
    ENDIF()

    # Add the library and all the dependencies

    IF (PARSE_DEFINES)
      ADD_DEFINITIONS(${PARSE_DEFINES})
    ENDIF()

    IF (NOT PARSE_CUDALIBRARY)
      ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ELSE()
      CUDA_ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
        ${PARSE_SOURCES})
    ENDIF()

    SET_PROPERTY(TARGET ${LIBRARY_NAME} APPEND PROPERTY
      LABELS ${PACKAGE_NAME})

    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

    TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${LINK_LIBS})

    # Add to the install target

    IF (NOT PARSE_TESTONLY AND NOT PARSE_NO_INSTALL_LIB_OR_HEADERS)

      SET_PROPERTY(GLOBAL PROPERTY ${PROJECT_NAME}_HAS_INSTALL_TARGETS ON)
      SET_PROPERTY(GLOBAL PROPERTY ${PACKAGE_NAME}_HAS_INSTALL_TARGETS ON)

      INSTALL(
        TARGETS ${LIBRARY_NAME}
        EXPORT ${PROJECT_NAME}
          RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
          LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
          ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )

      INSTALL(
        FILES ${PARSE_HEADERS}
        DESTINATION "${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}"
        COMPONENT ${PACKAGE_NAME}
        )

    ELSE()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping installation hooks for this library because 'TESTONLY' was passed in ...")
      ENDIF()

    ENDIF()

    # Append the new include dirs, library dirs, and libraries to this package's lists

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT PACKAGE_LIBRARY_DIRS)

    IF (NOT PARSE_TESTONLY)

      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
      IF (PARSE_IMPORTEDLIBS)
        PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES ${PARSE_IMPORTEDLIBS})
      ENDIF()
      PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES ${LIBRARY_NAME})

      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_INCLUDE_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARY_DIRS)
      REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARIES)

    ELSE()

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(STATUS "Skipping augmentation of package's lists of include"
          " directories and libraries because 'TESTONLY' was passed in ...")
      ENDIF()

      GLOBAL_SET(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        PRINT_VAR(${LIBRARY_NAME}_INCLUDE_DIRS)
      ENDIF()

    ENDIF()
  ENDIF() #if not in installation testing mode

  IF(${PROJECT_NAME}_ENABLE_INSTALLATION_TESTING)
    LIST(FIND Trilinos_INSTALLATION_PACKAGE_LIST ${PACKAGE_NAME} ${PACKAGE_NAME}_WAS_INSTALLED)
    IF(${${PACKAGE_NAME}_WAS_INSTALLED} EQUAL -1)
      MESSAGE(FATAL_ERROR
        "The package ${PACKAGE_NAME} was not installed with ${PROJECT_NAME}! Please disable package ${PACKAGE_NAME} or install it.")
    ENDIF()

    INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING BEFORE ${${PACKAGE_NAME}_INSTALLATION_INCLUDE_DIRS} ${${PACKAGE}_INSTALLATION_TPL_INCLUDE_DIRS})
    SET_PROPERTY(DIRECTORY APPEND PROPERTY PACKAGE_LIBRARY_DIRS ${${PACKAGE_NAME}_INSTALLATION_LIBRARY_DIRS})

    GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
    GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT PACKAGE_LIBRARY_DIRS)
    
    GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
    GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES    ${${PACKAGE_NAME}_INSTALLATION_LIBRARIES})

  ENDIF() #instalation testing mode

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()
  
ENDFUNCTION()
