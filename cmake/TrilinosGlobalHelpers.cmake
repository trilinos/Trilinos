
INCLUDE(Trilinos_Add_Option)
INCLUDE(TrilinosHelpers)
INCLUDE(Multiline_Set)
INCLUDE(Global_Null_Set)
INCLUDE(Remove_Global_Duplicates)
INCLUDE(Advanced_Option)

#
# 2008/10/06: rabartl:
#
# Below, we change the value of user cache values like
# Trilinos_ENABLE_${PACKAGE_NAME}, ${PACKAGE_NAME}_ENABLE_TESTS, and
# ${PACKAGE_NAME}_ENABLE_EXAMPLES by just setting them to regular variables
# instead of putting them back on the cache.  That means that they are used as
# global varibles but we don't want to disturb the cache since that would
# change the behavior for future invocations of cmake.  Because of this, these
# macros must all be called from the top-level Trilinos CMakeLists.txt file
# and macros must call macros as not to change the variable scope.
#
# I had to do it this way in order to be able to get the right behavior which
# is:
#
# 1) Override the value of these variables in all CMake processing
#
# 2) Avoid changing the user cache values because that would be confusing and
# would make it hard to change package enables/disable later without blowing
# away the cache
# 


#
# Macro that processes Trilinos_PACKAGES_AND_DIRS_AND_ENABLES into
# Trilinos_PACKAGES, Trilinos_PACKAGE_DIRS, Trilinos_NUM_PACKAGES,
# Trilinos_LAST_PACKAGE_IDX, and Trilinos_REVERSE_PACKAGES.
#
# This macro also sets up the standard package options along with
# default enables/disables.
#

MACRO(TRILINOS_PROCESS_PACKAGES_AND_DIRS_LISTS)

  ADVANCED_OPTION(Trilinos_REMOVE_DEFAULT_PACKAGE_DISABLES
    "Removes all default disables from the packages list.  Used for testing etc."
    OFF )

  # Separate out separate lists of package names and directoires
  
  LIST(LENGTH Trilinos_PACKAGES_AND_DIRS_AND_ENABLES Trilinos_NUM_PACKAGES_3)
  MATH(EXPR Trilinos_NUM_PACKAGES "${Trilinos_NUM_PACKAGES_3}/3")
  PRINT_VAR(Trilinos_NUM_PACKAGES)
  MATH(EXPR Trilinos_LAST_PACKAGE_IDX "${Trilinos_NUM_PACKAGES}-1")
  
  SET(Trilinos_PACKAGES)
  SET(Trilinos_PACKAGE_DIRS)
  FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
    MATH(EXPR PACKAGE_NAME_IDX "${PACKAGE_IDX}*3+0")
    MATH(EXPR PACKAGE_DIR_IDX "${PACKAGE_IDX}*3+1")
    MATH(EXPR PACKAGE_EANBLE_IDX "${PACKAGE_IDX}*3+2")
    LIST(GET Trilinos_PACKAGES_AND_DIRS_AND_ENABLES ${PACKAGE_NAME_IDX} PACKAGE)
    LIST(GET Trilinos_PACKAGES_AND_DIRS_AND_ENABLES ${PACKAGE_DIR_IDX} PACKAGE_DIR)
    LIST(GET Trilinos_PACKAGES_AND_DIRS_AND_ENABLES ${PACKAGE_EANBLE_IDX} PACKAGE_ENABLE)
    LIST(APPEND Trilinos_PACKAGES ${PACKAGE})
    LIST(APPEND Trilinos_PACKAGE_DIRS ${PACKAGE_DIR})
    IF (Trilinos_REMOVE_DEFAULT_PACKAGE_DISABLES AND "${PACKAGE_ENABLE}" STREQUAL "OFF")
      SET(PACKAGE_ENABLE "")
    ENDIF()
    TRILINOS_INSERT_STANDARD_PACKAGE_OPTIONS(${PACKAGE} "${PACKAGE_ENABLE}")
  ENDFOREACH()
  
  IF (Trilinos_VERBOSE_CONFIGURE)
    PRINT_VAR(Trilinos_PACKAGES)
    PRINT_VAR(Trilinos_PACKAGE_DIRS)
  ENDIF()
  
  # Create a reverse list for later use
  
  SET(Trilinos_REVERSE_PACKAGES ${Trilinos_PACKAGES})
  LIST(REVERSE Trilinos_REVERSE_PACKAGES)

ENDMACRO()


#
# Macro that reads all the package dependencies
#

MACRO(TRILINOS_READ_ALL_PACKAGE_DEPENDENCIES)
  
  FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
    LIST(GET Trilinos_PACKAGES ${PACKAGE_IDX} PACKAGE)
    LIST(GET Trilinos_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    TRILINOS_READ_PACKAGE_DEPENDENCIES(${PACKAGE} ${PACKAGE_DIR})
    TRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES(${PACKAGE})
  ENDFOREACH()
  
  ADVANCED_OPTION(Trilinos_DUMP_PACKAGE_DEPENDENCIES
    "Dump the package dependency information." OFF)
  
  IF (Trilinos_VERBOSE_CONFIGURE OR Trilinos_DUMP_PACKAGE_DEPENDENCIES)
    MESSAGE("")
    MESSAGE("Printing package dependencies ...")
    MESSAGE("")
    PRINT_VAR(Trilinos_PACKAGES)
    MESSAGE("")
    FOREACH(PACKAGE ${Trilinos_PACKAGES})
      TRILINOS_PRINT_PACKAGE_DEPENDENCIES(${PACKAGE})
      MESSAGE("")
    ENDFOREACH()
  ENDIF()

ENDMACRO()



#
# Macro that processes the list of TPLs
#

MACRO(TRILINOS_PROCESS_TPLS_LISTS)

  LIST(LENGTH Trilinos_TPLS_ENABLED Trilinos_NUM_TPLS_2)
  MATH(EXPR Trilinos_NUM_TPLS "${Trilinos_NUM_TPLS_2}/2")
  PRINT_VAR(Trilinos_NUM_TPLS)
  MATH(EXPR Trilinos_LAST_TPL_IDX "${Trilinos_NUM_TPLS}-1")
  
  SET(Trilinos_TPLS)
  FOREACH(TPL_IDX RANGE ${Trilinos_LAST_TPL_IDX})
    MATH(EXPR TPL_NAME_IDX "${TPL_IDX}*2")
    MATH(EXPR TPL_ENABLED_IDX "${TPL_IDX}*2+1")
    LIST(GET Trilinos_TPLS_ENABLED ${TPL_NAME_IDX} TPL)
    LIST(GET Trilinos_TPLS_ENABLED ${TPL_ENABLED_IDX} TPL_ENABLED)
    LIST(APPEND Trilinos_TPLS ${TPL})
    MULTILINE_SET(DOCSTR
      "Enable support for the TPL ${TPL} in all supported Trilinos packages."
      "  This can be set to 'ON', 'OFF', or left empty ''.  If Set to 'ON',"
      " then support for the TPL ${TPL} will be turned on by default in all packages"
      " that have optional support for this TPL.  However, the TPL will really"
      " only be enabled if there is at least one package having a dependency on the"
      " TPL.  If set to 'OFF', then all packages that have a required dependency on"
      " the TPL ${TPL} will be disabled and optional support for the TPL in all"
      " enabled packages will be turned off.  If left empty '', then other logic will"
      " be used to determine if the TPL ${TPL} should be enabled or not.  For example,"
      " if a package is enabled that has a required dependency on the TPL, then the"
      " TPL will be enabled and will be searched for.  If the TPL is enabled after"
      " all of this, then the user may need to set the location of the TPL in the"
      " variables TPL_${TPL}_INCLUDE_DIRS, TPL_${TPL}_LIBRARY_DIRS, and/or"
      " TPL_${TPL}_LIBRARIES.  For many of the TPLs, this information can be found"
      " out automatically by looking in standard locations."
      )
    SET(TPL_ENABLE_${TPL} ${TPL_ENABLED} CACHE STRING ${DOCSTR})
    # 2008/11/25: rabartl: Above, we use the prefix TPL_ instead of TRILINOS_ in order to
    # make it clear that external TPLs are different from Trilinos packages so users
    # don't get confused and think that Trilinos actually includes some TPL when it
    # does not!
  ENDFOREACH()
  
  IF (Trilinos_VERBOSE_CONFIGURE)
    PRINT_VAR(Trilinos_TPLS)
  ENDIF()
  
  # Create a reverse list for later use
  
  SET(Trilinos_REVERSE_TPLS ${Trilinos_TPLS})
  LIST(REVERSE Trilinos_REVERSE_TPLS)

ENDMACRO()


#
# Macro that gets the current list of enables packages
#
# Accesses the global varaibles:
#
#   Trilinos_PACKAGES
#   Trilinos_ENABLE_${PACKAGE}
#
# where ${PACKAGE} is every package in TRILINOS_PACKAGES_IN
#

MACRO(TRILINOS_GET_ENABLED_LIST
  LISTVAR ENABLED_PREFIX ENABLED_FLAG
  TRILINOS_ENABLED_LIST_OUT NUM_ENABLED_OUT
  )
  SET(${TRILINOS_ENABLED_LIST_OUT} "")
  SET(${NUM_ENABLED_OUT} 0)
  FOREACH(ENTITY ${${LISTVAR}})
    ASSERT_DEFINED(${ENABLED_PREFIX}_ENABLE_${ENTITY})
    IF ("${${ENABLED_PREFIX}_ENABLE_${ENTITY}}" STREQUAL "${ENABLED_FLAG}")
      SET(${TRILINOS_ENABLED_LIST_OUT} "${${TRILINOS_ENABLED_LIST_OUT}} ${ENTITY}")
      MATH(EXPR ${NUM_ENABLED_OUT} "${${NUM_ENABLED_OUT}}+1")
    ENDIF()
  ENDFOREACH()
ENDMACRO()


#
# Function that prints the current set of enabled/disabled packages
#

FUNCTION(TRILINOS_PRINT_ENABLED_PACKAGE_LIST DOCSTRING ENABLED_FLAG)
  TRILINOS_GET_ENABLED_LIST( Trilinos_PACKAGES Trilinos ${ENABLED_FLAG}
    Trilinos_ENABLED_PACKAGES NUM_ENABLED)
  MESSAGE("${DOCSTRING}: ${Trilinos_ENABLED_PACKAGES} ${NUM_ENABLED}")
ENDFUNCTION()


#
# Function that prints the current set of enabled/disabled TPLs
#

FUNCTION(TRILINOS_PRINT_ENABLED_TPL_LIST DOCSTRING ENABLED_FLAG)
  TRILINOS_GET_ENABLED_LIST( Trilinos_TPLS TPL ${ENABLED_FLAG}
    Trilinos_ENABLED_PACKAGES NUM_ENABLED)
  MESSAGE("${DOCSTRING}: ${Trilinos_ENABLED_PACKAGES} ${NUM_ENABLED}")
ENDFUNCTION()


#
# Function that sets a varaible to DECLARED-UNDEFINED
#

FUNCTION(DECLARE_UNDEFINED VAR)
  SET(${VAR} DECLARED-UNDEFINED PARENT_SCOPE)
ENDFUNCTION()


#
# Function that asserts that a package dependency variable is defined
# correctly
#

FUNCTION(ASSERT_DEFINED_PACKAGE_VAR PACKAGE_VAR PACKAGE_NAME)
  IF (${PACKAGE_VAR} STREQUAL DECLARED-UNDEFINED)
    MESSAGE(FATAL_ERROR
      "Error, the package variable ${PACKAGE_VAR} was not defined correctly for package ${PACKAGE_NAME}!"
      )
  ENDIF()
ENDFUNCTION()


#
# Macro that reads in package dependencies for a package and sets forward
# dependencies for packages already read in.
#
# Modifies the global varibles:
#
#   ${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES
#   ${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES
#

MACRO(TRILINOS_READ_PACKAGE_DEPENDENCIES PACKAGE_NAME PACKAGE_DIR)

  SET(${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES "")
  SET(${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES "")

  DECLARE_UNDEFINED(LIB_REQUIRED_DEP_PACKAGES)
  DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_PACKAGES)
  DECLARE_UNDEFINED(TEST_REQUIRED_DEP_PACKAGES)
  DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_PACKAGES)

  DECLARE_UNDEFINED(LIB_REQUIRED_DEP_TPLS "")
  DECLARE_UNDEFINED(LIB_OPTIONAL_DEP_TPLS "")
  DECLARE_UNDEFINED(TEST_REQUIRED_DEP_TPLS "")
  DECLARE_UNDEFINED(TEST_OPTIONAL_DEP_TPLS "")

  INCLUDE(packages/${PACKAGE_DIR}/cmake/Dependencies.cmake)

  ASSERT_DEFINED_PACKAGE_VAR(LIB_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(LIB_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(TEST_REQUIRED_DEP_PACKAGES ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(TEST_OPTIONAL_DEP_PACKAGES ${PACKAGE_NAME})

  ASSERT_DEFINED_PACKAGE_VAR(LIB_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(LIB_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(TEST_REQUIRED_DEP_TPLS ${PACKAGE_NAME})
  ASSERT_DEFINED_PACKAGE_VAR(TEST_OPTIONAL_DEP_TPLS ${PACKAGE_NAME})

  SET(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES ${LIB_REQUIRED_DEP_PACKAGES})
  SET(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES ${LIB_OPTIONAL_DEP_PACKAGES})
  SET(${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES ${TEST_REQUIRED_DEP_PACKAGES})
  SET(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES ${TEST_OPTIONAL_DEP_PACKAGES})

  SET(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS ${LIB_REQUIRED_DEP_TPLS})
  SET(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS ${LIB_OPTIONAL_DEP_TPLS})
  SET(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS ${TEST_REQUIRED_DEP_TPLS})
  SET(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS ${TEST_OPTIONAL_DEP_TPLS})

  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} LIB_REQUIRED_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} LIB_OPTIONAL_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} TEST_REQUIRED_DEP_PACKAGES)
  TRILINOS_APPEND_FORWARD_DEP_PACKAGES(${PACKAGE_NAME} TEST_OPTIONAL_DEP_PACKAGES)

ENDMACRO()


#
# Macro that prints out dependencies for a package
#
# Does not modify the global state.
#

MACRO(TRILINOS_PRINT_PACKAGE_DEPENDENCIES PACKAGE_NAME)

  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES)

  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_FORWARD_TEST_OPTIONAL_DEP_PACKAGES)

  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS)
  PRINT_NONEMPTY_VAR(${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS)

ENDMACRO()


#
# Macro that sets up standard user options for each package
#

MACRO(TRILINOS_INSERT_STANDARD_PACKAGE_OPTIONS PACKAGE_NAME PACKAGE_ENABLE)

  #MESSAGE("TRILINOS_INSERT_STANDARD_PACKAGE_OPTIONS: ${PACKAGE_NAME}")

  MULTILINE_SET(DOCSTR
    "Enable the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave"
    " empty to allow for other logic to decide.  If explicitly set to 'OFF'"
    " then all packages that have a required dependency on this package will"
    " be disabled and will override other explicit enables.   Look at the output"
    " from configure to see what packages are disabled in this case.  If set to 'ON',"
    " then the package ${PACKAGE_NAME} will be enabled if at all possible (i.e."
    " unless some other required package was explicitly disabled).  If left empty ''"
    " then other logic will be used to determine if the package will be enabled or"
    " disabled."
    )
  SET( Trilinos_ENABLE_${PACKAGE_NAME} "${PACKAGE_ENABLE}" CACHE STRING ${DOCSTR})

  MULTILINE_SET(DOCSTR
    "Build tests for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
     " to allow for other logic to decide.  This option only takes effect if the package"
     " is enabled after all of the enable/disable package logic has been applied."
     )
  SET( ${PACKAGE_NAME}_ENABLE_TESTS "" CACHE STRING ${DOCSTR})

  MULTILINE_SET(DOCSTR
    "Build examples for the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty ''"
     " to allow for other logic to decide.  This option only takes effect if the package"
     " is enabled after all of the enable/disable package logic has been applied."
     )
  SET( ${PACKAGE_NAME}_ENABLE_EXAMPLES "" CACHE STRING ${DOCSTR})

ENDMACRO()


#
# Macro that helps to set up forward package dependency lists
#

FUNCTION(TRILINOS_APPEND_FORWARD_DEP_PACKAGES PACKAGE_NAME LIST_TYPE)

  SET(DEP_PKG_LIST_NAME "${PACKAGE_NAME}_${LIST_TYPE}")

  #MESSAGE("DEP_PKG_LIST_NAME = ${DEP_PKG_LIST_NAME}")
  #MESSAGE("${DEP_PKG_LIST_NAME} = ${${DEP_PKG_LIST_NAME}}")

  FOREACH(DEP_PKG ${${DEP_PKG_LIST_NAME}})
    #MESSAGE("DEP_PKG = ${DEP_PKG}")
    SET(FWD_DEP_PKG_LIST_NAME "${DEP_PKG}_FORWARD_${LIST_TYPE}")
    #MESSAGE("FWD_DEP_PKG_LIST_NAME = ${FWD_DEP_PKG_LIST_NAME}")
    IF (NOT DEFINED ${FWD_DEP_PKG_LIST_NAME})
      MULTILINE_SET(ERRMSG
        "Error, the package '${DEP_PKG}' is listed as a dependency of the package"
        " '${PACKAGE_NAME}' in the list '${DEP_PKG_LIST_NAME}' but the package"
        " '${DEP_PKG}' is either not defined or is listed later in the package order."
        "  Check the spelling of '${DEP_PKG}' or see how it is listed in"
        " Trilinos_PACKAGES_AND_DIRS_AND_ENABLES in relationship to '${PACKAGE_NAME}'.")
      MESSAGE(FATAL_ERROR ${ERRMSG})
    ENDIF()
    SET(${FWD_DEP_PKG_LIST_NAME} ${${FWD_DEP_PKG_LIST_NAME}} ${PACKAGE_NAME} PARENT_SCOPE)
  ENDFOREACH()

ENDFUNCTION()


#
# Private helper macros
#

MACRO(TRILINOS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE PACKAGE_NAME OPTIONAL_DEP_PACKAGE)

  SET( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} "" CACHE STRING
    "Enable optional support for the package ${OPTIONAL_DEP_PACKAGE} in the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty to allow for other logic to decide" )

ENDMACRO()

MACRO(TRILINOS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE PACKAGE_NAME OPTIONAL_DEP_TPL)

  SET( ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} "" CACHE STRING
    "Enable optional support for the TPL ${OPTIONAL_DEP_TPL} in the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave empty to allow for other logic to decide" )

ENDMACRO()


#
# Macro that enables optional package interdependencies dependancies
#

MACRO(TRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nTRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    TRILINOS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
    TRILINOS_PRIVATE_ADD_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
    TRILINOS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
    TRILINOS_PRIVATE_ADD_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
  ENDFOREACH()

ENDMACRO()


#
# Private helper macro
#

MACRO(TRILINOS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME LIBRARY_DEP)

  #MESSAGE("TRILINOS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES ${FORWARD_DEP_PACKAGE_NAME} ${LIBRARY_DEP}")  

  IF ("${LIBRARY_DEP}" STREQUAL "TRUE")
    SET(DEP_TYPE_STR "library")
    ASSERT_DEFINED(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
    IF (Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME} OR Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME} STREQUAL "")
      MESSAGE(STATUS
        "Setting Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME}=OFF because ${FORWARD_DEP_PACKAGE_NAME} has a required ${DEP_TYPE_STR} dependence on disabled package ${PACKAGE_NAME}")
      SET(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME} OFF)
    ENDIF()
  ELSE()
    SET(DEP_TYPE_STR "test/example")
  ENDIF()

  ASSERT_DEFINED(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS)
  IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS OR ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS STREQUAL "")
    # Always disable the option but don't print the message if the package is not enabled
    IF (Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
      MESSAGE(STATUS
        "Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS=OFF because ${FORWARD_DEP_PACKAGE_NAME} has a required ${DEP_TYPE_STR} dependence on disabled package ${PACKAGE_NAME}")
    ENDIF()
    SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_TESTS OFF)
  ENDIF()

  ASSERT_DEFINED(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES)
  IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES OR ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES STREQUAL "")
    # Always disable the option but don't print the message if the package is not enabled
    IF (Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
      MESSAGE(STATUS
        "Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES=OFF because ${FORWARD_DEP_PACKAGE_NAME} has a required ${DEP_TYPE_STR} dependence on disabled package ${PACKAGE_NAME}")
    ENDIF()
    SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_EXAMPLES OFF)
  ENDIF()

ENDMACRO()


#
# Private helper macro
#

MACRO(TRILINOS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME)

  #MESSAGE("TRILINOS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES ${FORWARD_DEP_PACKAGE_NAME} ${PACKAGE_NAME}")  

  ASSERT_DEFINED(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME})
  IF (${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OR ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} STREQUAL "")
    # Always disable the conditional enable but only print the message if the package is enabled.
    IF (Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
      MESSAGE(STATUS
        "Setting ${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME}=OFF because ${FORWARD_DEP_PACKAGE_NAME} has an optional library dependence on disabled package ${PACKAGE_NAME}")
    ENDIF()
    SET(${FORWARD_DEP_PACKAGE_NAME}_ENABLE_${PACKAGE_NAME} OFF)
  ENDIF()

ENDMACRO()


#
# Function that disables all forward packages recurrsively
#

MACRO(TRILINOS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES PACKAGE_NAME)

  #MESSAGE("TRILINOS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES: ${PACKAGE_NAME}")

  IF ("${Trilinos_ENABLE_${PACKAGE}}" STREQUAL "OFF")

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      TRILINOS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME} TRUE)
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      TRILINOS_PRIVATE_DISABLE_OPTIONAL_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_TEST_REQUIRED_DEP_PACKAGES})
      TRILINOS_PRIVATE_DISABLE_REQUIRED_PACKAGE_ENABLES(${FWD_DEP_PKG} ${PACKAGE_NAME} FALSE)
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Private helper macros
#

MACRO(TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE PACKAGE_NAME OPTIONAL_DEP_PACKAGE)

  #MESSAGE("TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE}")

  IF("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}}" STREQUAL "")
    IF(Trilinos_ENABLE_${PACKAGE_NAME} AND Trilinos_ENABLE_${OPTIONAL_DEP_PACKAGE})
      MESSAGE(STATUS "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON since Trilinos_ENABLE_${PACKAGE_NAME}=ON AND Trilinos_ENABLE_${OPTIONAL_DEP_PACKAGE}=ON")
      SET(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE} ON)
    ENDIF()
  ENDIF()

  STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  STRING(TOUPPER ${OPTIONAL_DEP_PACKAGE} OPTIONAL_DEP_PACKAGE_UPPER)
  SET(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_PACKAGE_UPPER})

  IF(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_PACKAGE})
    SET(${MACRO_DEFINE_NAME} ON)
  ELSE()
    SET(${MACRO_DEFINE_NAME} OFF)
  ENDIF()

ENDMACRO()


MACRO(TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE PACKAGE_NAME OPTIONAL_DEP_TPL)

  #MESSAGE("TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE: ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL}")

  IF("${${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}}" STREQUAL "")
    IF(Trilinos_ENABLE_${PACKAGE_NAME} AND TPL_ENABLE_${OPTIONAL_DEP_TPL})
      MESSAGE(STATUS "Setting ${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL}=ON since TPL_ENABLE_${OPTIONAL_DEP_TPL}=ON")
      SET(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL} ON)
    ENDIF()
  ENDIF()

  STRING(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UPPER)
  STRING(TOUPPER ${OPTIONAL_DEP_TPL} OPTIONAL_DEP_TPL_UPPER)
  SET(MACRO_DEFINE_NAME HAVE_${PACKAGE_NAME_UPPER}_${OPTIONAL_DEP_TPL_UPPER})

  IF(${PACKAGE_NAME}_ENABLE_${OPTIONAL_DEP_TPL})
    SET(${MACRO_DEFINE_NAME} ON)
  ELSE()
    SET(${MACRO_DEFINE_NAME} OFF)
  ENDIF()

ENDMACRO()


#
# Macro that post-processes optional dependancies after all other
# dependencies have been worked out
#

MACRO(TRILINOS_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nTRILINOS_ADD_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
    TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_PACKAGE ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
    TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_PACKAGE_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_PACKAGE} )
  ENDFOREACH()

ENDMACRO()


#
# Macro that post-processes optional package TPL based on if the TPL
# has been enabled or not
#

MACRO(TRILINOS_POSTPROCESS_OPTIONAL_TPL_ENABLES PACKAGE_NAME)

  #MESSAGE("\nTRILINOS_ADD_OPTIONAL_TPL_ENABLES: ${PACKAGE_NAME}")

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_TPLS})
    TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
  ENDFOREACH()

  FOREACH(OPTIONAL_DEP_TPL ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_TPLS})
    TRILINOS_PRIVATE_POSTPROCESS_OPTIONAL_TPL_ENABLE(
      ${PACKAGE_NAME} ${OPTIONAL_DEP_TPL} )
  ENDFOREACH()

ENDMACRO()


#
# Set an individual pacakge variable based on the global value
#

MACRO(TRILINOS_POSTPROCESS_STANDARD_PACKAGE_VARIABLE TRILINOS_VAR PACKAGE_VAR)

  IF (Trilinos_VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRILINOS_POSTPROCESS_STANDARD_PACKAGE_VARIABLE:")
    MESSAGE(STATUS "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
    MESSAGE(STATUS "${TRILINOS_VAR} = ${${TRILINOS_VAR}}")
  ENDIF()

  IF (${PACKAGE_VAR} STREQUAL "")
    IF (${TRILINOS_VAR} STREQUAL "ON")
      MESSAGE(STATUS "Setting ${PACKAGE_VAR}=ON")
      SET(${PACKAGE_VAR} ON)
    ELSEIF (TRILINOS_VAR STREQUAL "OFF")
      MESSAGE(STATUS "Setting ${PACKAGE_VAR}=OFF")
      SET(${PACKAGE_VAR} OFF)
    ELSE()
      #MESSAGE(STATUS "ELSE")
      # Otherwise, we will leave it up the the individual package
      # to decide?
    ENDIF()
  ELSE()
    #MESSAGE(STATUS "PACKAGE_VAR NOT DEFAULT")
  ENDIF()

  IF (Trilinos_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  ENDIF()

ENDMACRO()


#
# Macro used to set Trilinos_ENABLE_${PACKAGE_NAME} based on
# Trilinos_ENABLE_ALL_PACKAGES
#

MACRO(TRILINOS_APPLY_ALL_PACKAGE_ENABLES PACKAGE_NAME)

  TRILINOS_POSTPROCESS_STANDARD_PACKAGE_VARIABLE(
    Trilinos_ENABLE_ALL_PACKAGES Trilinos_ENABLE_${PACKAGE_NAME} )

ENDMACRO()


#
# Macro used to set ${PACKAGE)_ENABLE_TESTS and ${PACKAGE)_ENABLE_EXAMPLES
# based on Trilinos_ENABLE_ALL_PACKAGES
#

MACRO(TRILINOS_APPLY_TEST_EXAMPLE_ENABLES PACKAGE_NAME)

  IF (Trilinos_ENABLE_${PACKAGE_NAME})

    TRILINOS_POSTPROCESS_STANDARD_PACKAGE_VARIABLE(
      Trilinos_ENABLE_TESTS ${PACKAGE_NAME}_ENABLE_TESTS )

    TRILINOS_POSTPROCESS_STANDARD_PACKAGE_VARIABLE(
      Trilinos_ENABLE_EXAMPLES ${PACKAGE_NAME}_ENABLE_EXAMPLES )

  ENDIF()

ENDMACRO()


#
# Private helper macro
#

MACRO(TRILINOS_PRIVATE_ENABLE_FORWARD_PACKAGE FORWARD_DEP_PACKAGE_NAME PACKAGE_NAME)
  # Enable the forward package if it is not already set to ON or OFF
  ASSERT_DEFINED(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
  IF(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME} STREQUAL "")
    MESSAGE(STATUS "Setting Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME}=ON because Trilinos_ENABLE_${PACKAGE_NAME}=ON")
    ASSERT_DEFINED(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME})
    SET(Trilinos_ENABLE_${FORWARD_DEP_PACKAGE_NAME} ON)
  ENDIF()
ENDMACRO()


#
# Macro used to set Trilinos_ENABLE_${FWD_PACKAGE_NAME)=ON for all optional
# and required forward dependencies of the package ${PACKAGE_NAME}
#

MACRO(TRILINOS_ENABLE_FORWARD_PACKAGE_ENABLES PACKAGE_NAME)

  #MESSAGE("\nTRILINOS_ENABLE_FORWARD_PACKAGE_ENABLES ${PACKAGE_NAME}")
  #MESSAGE(STATUS "Trilinos_ENABLE_${PACKAGE_NAME}=${Trilinos_ENABLE_${PACKAGE_NAME}}")

  # Enable the forward packages if this package is enabled
  ASSERT_DEFINED(Trilinos_ENABLE_${PACKAGE_NAME})
  IF (Trilinos_ENABLE_${PACKAGE_NAME})

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_REQUIRED_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

    FOREACH(FWD_DEP_PKG ${${PACKAGE_NAME}_FORWARD_LIB_OPTIONAL_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_FORWARD_PACKAGE(${FWD_DEP_PKG} ${PACKAGE_NAME})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Private helper macros
#

MACRO(TRILINOS_PRIVATE_ENABLE_DEP_PACKAGE PACKAGE_NAME DEP_PACKAGE_NAME)
  ASSERT_DEFINED(Trilinos_ENABLE_${DEP_PACKAGE_NAME})
  IF(Trilinos_ENABLE_${DEP_PACKAGE_NAME} STREQUAL "")
    MESSAGE(STATUS "Setting Trilinos_ENABLE_${DEP_PACKAGE_NAME}=ON because Trilinos_ENABLE_${PACKAGE_NAME}=ON")
    ASSERT_DEFINED(Trilinos_ENABLE_${DEP_PACKAGE_NAME})
    SET(Trilinos_ENABLE_${DEP_PACKAGE_NAME} ON)
  ENDIF()
ENDMACRO()

MACRO(TRILINOS_PRIVATE_ENABLE_DEP_TPL PACKAGE_NAME DEP_TPL_NAME)
  ASSERT_DEFINED(TPL_ENABLE_${DEP_TPL_NAME})
  IF(TPL_ENABLE_${DEP_TPL_NAME} STREQUAL "")
    MESSAGE(STATUS "Setting TPL_ENABLE_${DEP_TPL_NAME}=ON because it is required by the enabled package ${PACKAGE_NAME}")
    ASSERT_DEFINED(TPL_ENABLE_${DEP_TPL_NAME})
    SET(TPL_ENABLE_${DEP_TPL_NAME} ON)
  ENDIF()
ENDMACRO()


#
# Macro that sets the optional packages for given package
#

MACRO(TRILINOS_ENABLE_OPTIONAL_PACKAGES PACKAGE_NAME)

  #MESSAGE("TRILINOS_ENABLE_OPTIONAL_PACKAGE_ENABLES: ${PACKAGE_NAME}")
  #MESSAGE(STATUS "Trilinos_ENABLE_${PACKAGE_NAME}=${Trilinos_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(Trilinos_ENABLE_${PACKAGE_NAME})

  IF (Trilinos_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_OPTIONAL_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG})
    ENDFOREACH()

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_TEST_OPTIONAL_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that sets the required packages for given package
#

MACRO(TRILINOS_ENABLE_REQUIRED_PACKAGES PACKAGE_NAME)

  #MESSAGE("TRILINOS_ENABLE_REQUIRED_PACKAGE_ENABLES: ${PACKAGE_NAME}")
  #MESSAGE(STATUS "Trilinos_ENABLE_${PACKAGE_NAME}=${Trilinos_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(Trilinos_ENABLE_${PACKAGE_NAME})

  IF (Trilinos_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG})
    ENDFOREACH()

    FOREACH(DEP_PKG ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_PACKAGES})
      TRILINOS_PRIVATE_ENABLE_DEP_PACKAGE(${PACKAGE_NAME} ${DEP_PKG})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Macro that sets the required TPLs for given package
#

MACRO(TRILINOS_ENABLE_REQUIRED_TPLS PACKAGE_NAME)

  #MESSAGE("TRILINOS_ENABLE_REQUIRED_TPL_ENABLES: ${PACKAGE_NAME}")
  #MESSAGE(STATUS "Trilinos_ENABLE_${PACKAGE_NAME}=${Trilinos_ENABLE_${PACKAGE_NAME}}")

  ASSERT_DEFINED(Trilinos_ENABLE_${PACKAGE_NAME})

  IF (Trilinos_ENABLE_${PACKAGE_NAME})

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_LIB_REQUIRED_DEP_TPLS})
      TRILINOS_PRIVATE_ENABLE_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

    FOREACH(DEP_TPL ${${PACKAGE_NAME}_TEST_REQUIRED_DEP_TPLS})
      TRILINOS_PRIVATE_ENABLE_DEP_TPL(${PACKAGE_NAME} ${DEP_TPL})
    ENDFOREACH()

  ENDIF()

ENDMACRO()


#
# Private helper stuff
#


FUNCTION(TRILINOS_WRITE_DEPS_TO_XML_FILE PACKAGE_NAME LIST_TYPE)

  SET(DEPS_VAR ${PACKAGE_NAME}_${LIST_TYPE})
  ASSERT_DEFINED(DEPS_VAR)
  SET(DEPS ${${DEPS_VAR}})

  #PRINT_VAR(PACKAGE_NAME)
  #PRINT_VAR(DEPS)

  IF (NOT DEPS)

    FILE(APPEND ${Trilinos_DEPS_XML_OUTPUT_FILE}
      "    <${LIST_TYPE}/>\n")
    
  ELSE()

    SET(VALUE_STR "")

    FOREACH(DEP ${DEPS})

      IF(VALUE_STR)
        SET(VALUE_STR "${VALUE_STR},")
      ENDIF()

      SET(VALUE_STR "${VALUE_STR}${DEP}")

    ENDFOREACH()

    FILE(APPEND ${Trilinos_DEPS_XML_OUTPUT_FILE}
      "    <${LIST_TYPE} value=\"${VALUE_STR}\"/>\n")

  ENDIF()

ENDFUNCTION()


#
# Function that writes the dependency information for Trilinos into
# an XML file for other tools to use.
#

FUNCTION(TRILINOS_DUMP_DEPS_XML_FILE)

  FILE(WRITE ${Trilinos_DEPS_XML_OUTPUT_FILE}
    "<PackageDependencies>\n" )

  FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
    LIST(GET Trilinos_PACKAGES ${PACKAGE_IDX} PACKAGE)
    LIST(GET Trilinos_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    
    FILE(APPEND ${Trilinos_DEPS_XML_OUTPUT_FILE}
      "  <Package name=\"${PACKAGE}\" dir=\"${PACKAGE_DIR}\">\n")

    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} LIB_REQUIRED_DEP_PACKAGES)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} LIB_OPTIONAL_DEP_PACKAGES)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} TEST_REQUIRED_DEP_PACKAGES)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} TEST_OPTIONAL_DEP_PACKAGES)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} LIB_REQUIRED_DEP_TPLS)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} LIB_OPTIONAL_DEP_TPLS)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} TEST_REQUIRED_DEP_TPLS)
    TRILINOS_WRITE_DEPS_TO_XML_FILE(${PACKAGE} TEST_OPTIONAL_DEP_TPLS)

    FILE(APPEND ${Trilinos_DEPS_XML_OUTPUT_FILE}
      "  </Package>\n")

  ENDFOREACH()

  FILE(APPEND ${Trilinos_DEPS_XML_OUTPUT_FILE}
    "</PackageDependencies>\n" )
  

ENDFUNCTION()


#
# Macro that adjusts all of the package enables from what the user input
# to the final set that will be used to enable packages
#

MACRO(TRILINOS_ADJUST_PACKAGE_ENABLES)

  MESSAGE("")
  MESSAGE("Disabling all packages that have a required dependency on explicitly disabled TPLs ...")
  MESSAGE("")
  # ToDO: Implement This!

  MESSAGE("")
  MESSAGE("Disabling forward required packages and optional intra-package support that have a dependancy on explicitly disabled packages ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_DISABLE_FORWARD_REQUIRED_DEP_PACKAGES(${PACKAGE})
  ENDFOREACH()
  
  IF (Trilinos_ENABLE_ALL_PACKAGES)
    MESSAGE("")
    MESSAGE("Enabling all packages that are not currently disabled ...")
    MESSAGE("")
    FOREACH(PACKAGE ${Trilinos_PACKAGES})
      TRILINOS_APPLY_ALL_PACKAGE_ENABLES(${PACKAGE})
    ENDFOREACH()
  ENDIF()
  
  IF (Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES)
    MESSAGE("")
    MESSAGE("Enabling all forward dependent packages ...")
    MESSAGE("")
    FOREACH(PACKAGE ${Trilinos_PACKAGES})
      TRILINOS_ENABLE_FORWARD_PACKAGE_ENABLES(${PACKAGE})
    ENDFOREACH()
    SET(Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES ON)
  ENDIF()
  
  IF (Trilinos_ENABLE_TESTS OR Trilinos_ENABLE_EXAMPLES)
    MESSAGE("")
    MESSAGE("Enabling all tests and examples that have not been explicitly disabled ...")
    MESSAGE("")
    FOREACH(PACKAGE ${Trilinos_PACKAGES})
      TRILINOS_APPLY_TEST_EXAMPLE_ENABLES(${PACKAGE})
    ENDFOREACH()
  ENDIF()
  # NOTE: Above, we enable tests and examples here, before the remaining required
  # packages so that we don't enable tests that don't need to be enabled based
  # on the use of the option Trilinos_ENABLE_ALL_FORWARD_DEP_PACKAGES.
  
  IF (Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES)
    MESSAGE("")
    MESSAGE("Enabling all optional packages for current set of enabled packages ...")
    MESSAGE("")
    FOREACH(PACKAGE ${Trilinos_REVERSE_PACKAGES})
      TRILINOS_ENABLE_OPTIONAL_PACKAGES(${PACKAGE})
    ENDFOREACH()
  ENDIF()
  # NOTE: Above, we have to loop through the packages backward to enable all the
  # packages that feed into these packages.
  # NOTE Above, we don't have to enable the required packages because that will
  # come next
  
  MESSAGE("")
  MESSAGE("Enabling all remaining required packages for the current set of enabled packages ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_REVERSE_PACKAGES})
    TRILINOS_ENABLE_REQUIRED_PACKAGES(${PACKAGE})
  ENDFOREACH()
  
  MESSAGE("")
  MESSAGE("Enabling all optional intra-package enables that are not explicitly disabled if both sets of packages are enabled ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_POSTPROCESS_OPTIONAL_PACKAGE_ENABLES(${PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Enabling all remaining required TPLs for current set of enabled packages ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_ENABLE_REQUIRED_TPLS(${PACKAGE})
  ENDFOREACH()

  MESSAGE("")
  MESSAGE("Enabling all optional package TPL support for currently enabled TPLs ...")
  MESSAGE("")
  FOREACH(PACKAGE ${Trilinos_PACKAGES})
    TRILINOS_POSTPROCESS_OPTIONAL_TPL_ENABLES(${PACKAGE})
  ENDFOREACH()

ENDMACRO()



#
# Macro that gathers information from enabled TPLs
#

MACRO(TRILINOS_PROCESS_ENABLED_TPLS)
  FOREACH(TPL ${Trilinos_TPLS})
    IF (TPL_ENABLE_${TPL})
      MESSAGE(STATUS "Processing enabled TPL: ${TPL}")
      INCLUDE(TPLs/FindTPL${TPL})
    ENDIF()
  ENDFOREACH()
ENDMACRO()



#
# Macro that defines Trilinos testing support
#

MACRO(TRILINOS_SETUP_TESTING_SUPPORT)

  INCLUDE(CTest)
  
  IF (WIN32 AND NOT CYGWIN)
    SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT OFF)
  ELSE()
    SET(Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT ON)
  ENDIF()
  
  ADVANCED_OPTION(Trilinos_ENABLE_NATIVE_TEST_HARNESS
    "Enable the native Trilinos perl-based test harness."
    ${Trilinos_ENABLE_NATIVE_TEST_HARNESS_DEFAULT} )
  
  IF (Trilinos_ENABLE_NATIVE_TEST_HARNESS)
  
    ADD_CUSTOM_TARGET(
      runtests-serial
       ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
      --trilinos-dir=${TRILINOS_HOME_DIR}
      --comm=serial
      --build-dir=${TRILINOS_BUILD_DIR}
      --category=${TRILINOS_TEST_CATEGORY}
      --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
      )
  
    IF (TPL_ENABLE_MPI)
    
      ADD_CUSTOM_TARGET(
        runtests-mpi
         ${PERL_EXECUTABLE} ${TRILINOS_HOME_DIR}/commonTools/test/utilities/runtests
        --trilinos-dir=${TRILINOS_HOME_DIR}
        --comm=mpi
        --mpi-go="${TRILINOS_MPI_GO}"
        --max-proc=${MPIEXEC_MAX_NUMPROCS}
        --build-dir=${TRILINOS_BUILD_DIR}
        --category=${TRILINOS_TEST_CATEGORY}
        --output-dir=${TRILINOS_BUILD_DIR}/runtests-results
        )
  
    ENDIF()
  
  ENDIF()
  
  IF (WIN32)
    SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT OFF)
  ELSE()
    SET(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT ON)
  ENDIF()
  
  # 2008/10/17: rabartl: Above, I can not turn these tests on by default
  # with cygwin because the custom script target is not working for some
  # reason.
  
  ADVANCED_OPTION(Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS
    "Enable dependency unit tests."
    ${Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS_DEFAULT}
    )
  
  IF (Trilinos_ENABLE_DEPENCENCY_UNIT_TESTS)
    ADD_SUBDIRECTORY(cmake/DependencyUnitTests)
    ADD_SUBDIRECTORY(cmake/python/UnitTests)
  ENDIF()

  CONFIGURE_FILE(
    ${Trilinos_SOURCE_DIR}/cmake/ctest/CTestCustom.ctest.in
    ${Trilinos_BINARY_DIR}/CTestCustom.ctest
    )

ENDMACRO()



#
# Macro that defines Trilinos packaging options:
#

MACRO(TRILINOS_DEFINE_PACKAGING)

  SET(CPACK_PACKAGE_DESCRIPTION "Trilinos provides algorithms and technologies for the solution of large-scale, complex multi-physics engineering and scientific problems.")
  SET(CPACK_PACKAGE_FILE_NAME "trilinos-setup-${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_INSTALL_DIRECTORY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_REGISTRY_KEY "Trilinos ${Trilinos_VERSION}")
  SET(CPACK_PACKAGE_NAME "trilinos")
  SET(CPACK_PACKAGE_VENDOR "Sandia National Laboratories")
  SET(CPACK_PACKAGE_VERSION "${Trilinos_VERSION}")
  SET(CPACK_SOURCE_GENERATOR "TGZ;TBZ2")
  SET(CPACK_SOURCE_FILE_NAME "trilinos-source-${Trilinos_VERSION}")
  
  IF(WIN32)
    SET(CPACK_GENERATOR "NSIS")
    SET(CPACK_NSIS_MODIFY_PATH ON)
  ENDIF()
  
  INCLUDE(CPack)

ENDMACRO()



#
# Macro that does the final set of package configurations
#

MACRO(TRILINOS_CONFIGURE_ENABLED_PACKAGES)

  GLOBAL_NULL_SET(Trilinos_INCLUDE_DIRS)
  GLOBAL_NULL_SET(Trilinos_LIBRARY_DIRS)
  GLOBAL_NULL_SET(Trilinos_LIBRARIES)
  
  FOREACH(PACKAGE_IDX RANGE ${Trilinos_LAST_PACKAGE_IDX})
    LIST(GET Trilinos_PACKAGES ${PACKAGE_IDX} PACKAGE)
    LIST(GET Trilinos_PACKAGE_DIRS ${PACKAGE_IDX} PACKAGE_DIR)
    IF(Trilinos_ENABLE_${PACKAGE})
      SET(PACKAGE_NAME_GLOBAL ${PACKAGE}) # For consistency checking
      IF (NOT EXISTS ${PROJECT_SOURCE_DIR}/packages/${PACKAGE_DIR}/CMakeLists.txt)
        MESSAGE(FATAL_ERROR "Error, the file ${PACKAGE_DIR}/CMakeLists.txt does not exist!")
      ENDIF()
      ADD_SUBDIRECTORY(packages/${PACKAGE_DIR})
      LIST(APPEND Trilinos_INCLUDE_DIRS ${${PACKAGE}_INCLUDE_DIRS})
      LIST(APPEND Trilinos_LIBRARY_DIRS ${${PACKAGE}_LIBRARY_DIRS})
      LIST(APPEND Trilinos_LIBRARIES ${${PACKAGE}_LIBRARIES})
    ENDIF()
  ENDFOREACH()
  
  REMOVE_GLOBAL_DUPLICATES(Trilinos_INCLUDE_DIRS)
  REMOVE_GLOBAL_DUPLICATES(Trilinos_LIBRARY_DIRS)
  REMOVE_GLOBAL_DUPLICATES(Trilinos_LIBRARIES)
  
ENDMACRO()
