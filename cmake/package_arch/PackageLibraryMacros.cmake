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

MACRO(PACKAGE_CONFIGURE_FILE PACKAGE_NAME_CONFIG_FILE)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  ENDIF()

  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

ENDMACRO()


#
# Macro used to add a package library
#

FUNCTION(PACKAGE_ADD_LIBRARY LIBRARY_NAME)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nPACKAGE_ADD_LIBRARY: ${LIBRARY_NAME}")
  ENDIF()

  PARSE_ARGUMENTS(
    PARSE #prefix
    "HEADERS;NOINSTALLHEADERS;SOURCES;DEPLIBS;IMPORTEDLIBS" # Lists
    "TESTONLY" #Options
    ${ARGN} # Remaining arguments passed in
    )

  # Add the link directory for this library.

  LINK_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

  # NOTE: Above , this link path not really used here for anything.
  # Instead it is just added to the other set link library directories
  # that are already set.  These link directories are then extracted
  # and stored into stored in ${PACKAGE_NAME}_LIBRARY_DIRS.

  # Add whatever include directories have been defined so far

  INCLUDE_DIRECTORIES(AFTER ${${PACKAGE_NAME}_INCLUDE_DIRS})

  # Add whatever link directories have been added so far

  LINK_DIRECTORIES(${${PACKAGE_NAME}_LIBRARY_DIRS})
  LINK_DIRECTORIES(${${PACKAGE_NAME}_TEST_LIBRARY_DIRS})

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

  ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
    ${PARSE_SOURCES})

  SET_PROPERTY(TARGET ${LIBRARY_NAME} APPEND PROPERTY
    LABELS ${PACKAGE_NAME})

  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})
  
  TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${LINK_LIBS})

  # Add to the install target

  IF (NOT PARSE_TESTONLY)

    INSTALL(
      TARGETS ${LIBRARY_NAME}
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
  GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)

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

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

ENDFUNCTION()
