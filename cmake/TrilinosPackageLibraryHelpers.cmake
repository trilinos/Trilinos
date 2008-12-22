INCLUDE(Parse_Variable_Arguments)
INCLUDE(Global_Set)
INCLUDE(Append_Set)
INCLUDE(Append_Global_Set)
INCLUDE(Prepend_Global_Set)
INCLUDE(Remove_Global_Duplicates)
INCLUDE(TrilinosGeneralHelpers)


#
# Macro that configures the package's main config.h file
#

MACRO(TRILINOS_PACKAGE_CONFIGURE_FILE PACKAGE_NAME_CONFIG_FILE)

  IF (Trilinos_VERBOSE_CONFIGURE)
    MESSAGE("\nTRILINOS_PACKAGE_CONFIGURE_FILE: ${PACKAGE_NAME_CONFIG_FILE}")
  ENDIF()

  CONFIGURE_FILE(
    ${PACKAGE_SOURCE_DIR}/cmake/${PACKAGE_NAME_CONFIG_FILE}.in 
    ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    )

  INSTALL(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME_CONFIG_FILE}
    DESTINATION ${TRILINOS_INSTALL_LIB_INCLUDE_DIR}
    COMPONENT ${PACKAGE_NAME}
    )

ENDMACRO()


#
# Macro used to add a package library
#

FUNCTION(TRILINOS_PACKAGE_ADD_LIBRARY LIBRARY_NAME)

  IF (Trilinos_VERBOSE_CONFIGURE)
    MESSAGE("\nTRILINOS_PACKAGE_ADD_LIBRARY: ${LIBRARY_NAME}")
  ENDIF()

  PARSE_ARGUMENTS(
    PARSE #prefix
    "HEADERS;NOINSTALLHEADERS;SOURCES;DEPLIBS" # Lists
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

  IF (PARSE_DEPLIBS AND Trilinos_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "DEPLIBS = ${PARSE_DEPLIBS}")
  ENDIF()

  # Add the library and all the dependencies

  ADD_LIBRARY(${LIBRARY_NAME} ${PARSE_HEADERS} ${PARSE_NOINSTALLHEADERS}
    ${PARSE_SOURCES})

  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})

  SET(LINK_LIBS)

  IF (PARSE_DEPLIBS)
    APPEND_SET(LINK_LIBS ${PARSE_DEPLIBS})
  ENDIF()

  # We only want to link to the dependent package and TPL libraries when we need
  # to.  We only need to link to these dependent libraries when this is the first
  # library being created for this package or if this library does not depend
  # on other libraries created for this package.  Otherwise, we don't need to
  # add the include directories or link libraries because a dependent lib
  # specified in PARSE_DEP_LIBS already has everything that we need.  We also
  # Need to make special considerations for test libraries since things
  # need to be handled a little bit differently (but not much).

  SET(ADD_DEP_PACKAGE_AND_TPL_LIBS TRUE)

  IF (PARSE_DEPLIBS)
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
  
    # Add the dependent package libraries (if we have not done so yet for this package)
    PACKAGE_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} LIB PACKAGES ALL_DEP_PACKAGES)
    PACKAGE_SORT_AND_APPEND_PATHS_LIBS("${Trilinos_REVERSE_PACKAGES}" "${ALL_DEP_PACKAGES}"
      "" LINK_LIBS)

    # Add the TPL libraries (if we have not done so yet for this package)
    PACKAGE_GATHER_ENABLED_ITEMS(${PACKAGE_NAME} LIB TPLS ALL_TPLS)
    PACKAGE_SORT_AND_APPEND_PATHS_LIBS("${Trilinos_REVERSE_TPLS}" "${ALL_TPLS}"
      TPL_ LINK_LIBS)

  ENDIF()

  IF (Trilinos_VERBOSE_CONFIGURE)
    PRINT_VAR(LINK_LIBS)
  ENDIF()
  
  TARGET_LINK_LIBRARIES(${LIBRARY_NAME} ${LINK_LIBS})

  # Add to the install target

  IF (NOT PARSE_TESTONLY)

    INSTALL(
      TARGETS ${LIBRARY_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
      COMPONENT ${PACKAGE_NAME}
      )
    
    INSTALL(
      FILES ${PARSE_HEADERS}
      DESTINATION ${TRILINOS_INSTALL_INCLUDE_DIR}
      COMPONENT ${PACKAGE_NAME}
      )

  ELSE()

    IF (Trilinos_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Skipping installation hooks for this library because 'TESTONLY' was passed in ...")
    ENDIF()

  ENDIF()

  # Append the new include dirs, library dirs, and libraries to this package's lists

  GET_DIRECTORY_PROPERTY(INCLUDE_DIRS_CURRENT INCLUDE_DIRECTORIES)
  GET_DIRECTORY_PROPERTY(LIBRARY_DIRS_CURRENT LINK_DIRECTORIES)

  IF (NOT PARSE_TESTONLY)

    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS ${LIBRARY_DIRS_CURRENT})
    PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES ${LIBRARY_NAME})
  
    REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_INCLUDE_DIRS)
    REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARY_DIRS)
    REMOVE_GLOBAL_DUPLICATES(${PACKAGE_NAME}_LIBRARIES)

  ELSE()

    IF (Trilinos_VERBOSE_CONFIGURE)
      MESSAGE(STATUS "Skipping augmentation of package's lists of include directories and libraries because 'TESTONLY' was passed in ...")
    ENDIF()

    GLOBAL_SET(${LIBRARY_NAME}_INCLUDE_DIRS ${INCLUDE_DIRS_CURRENT})

    IF (Trilinos_VERBOSE_CONFIGURE)
      PRINT_VAR(${LIBRARY_NAME}_INCLUDE_DIRS)
    ENDIF()

  ENDIF()

  IF (Trilinos_VERBOSE_CONFIGURE)
    PRINT_VAR(${PACKAGE_NAME}_INCLUDE_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARY_DIRS)
    PRINT_VAR(${PACKAGE_NAME}_LIBRARIES)
  ENDIF()

  # 2008/11/21: rabartl: ToDo: Above: Generalize this for test-only
  # libraries!

ENDFUNCTION()


MACRO(TRILINOS_PACKAGE_EXPORT_DEPENDENCY_VARIABLES)
  # 2008/11/21: rabartl: ToDo: Get rid if this macro all together!
ENDMACRO()
