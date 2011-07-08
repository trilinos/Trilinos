
INCLUDE(SetCacheOnOffEmpty)
INCLUDE(MultilineSet)
INCLUDE(AdvancedOption)
INCLUDE(PackageListHelpers)


#
# Helper functions
#

MACRO(ALLOW_MISSING_EXTERNAL_PACKAGES)
  FOREACH(PACKAGE ${ARGN})
    SET(${PACKAGE}_ALLOW_MISSING_EXTERNAL_PACKAGE TRUE)
  ENDFOREACH()
ENDMACRO()



#
# Below, we change the value of user cache values like
# ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME},
# ${PACKAGE_NAME}_ENABLE_TESTS, and ${PACKAGE_NAME}_ENABLE_EXAMPLES by
# just setting them to regular variables that live in the top scope
# instead of putting them back in the cache.  That means that they are
# used as global variables but we don't want to disturb the cache
# since that would change the behavior for future invocations of cmake
# (which is very confusing).  Because of this, these macros must all
# be called from the top-level ${PROJECT_NAME} CMakeLists.txt file and
# macros must call macros as not to change the variable scope.
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
# Macro that sets up standard user options for each package
#

MACRO(PACKAGE_ARCH_INSERT_STANDARD_PACKAGE_OPTIONS  PACKAGE_NAME  PACKAGE_CLASSIFICATION)

  #MESSAGE("PACKAGE_ARCH_INSERT_STANDARD_PACKAGE_OPTIONS: ${PACKAGE_NAME} ${PACKAGE_CLASSIFICATION}")

  IF (${PACKAGE_CLASSIFICATION} STREQUAL PS OR ${PACKAGE_CLASSIFICATION} STREQUAL SS) 
    #MESSAGE("PS or SS")
    SET(PACKAGE_ENABLE "")
  ELSEIF (${PACKAGE_CLASSIFICATION} STREQUAL EX)
    #MESSAGE("EX")
    SET(PACKAGE_ENABLE OFF)
  ELSE()
    MESSAGE(FATAL_ERROR "Error the package classification '${PACKAGE_CLASSIFICATION}'"
      " for the package ${PACKAGE_NAME} is not a valid classification." )
  ENDIF()
  #PRINT_VAR(PACKAGE_ENABLE)

  IF (NOT ${PACKAGE_NAME}_CLASSIFICATION) # Allow testing override
    SET(${PACKAGE_NAME}_CLASSIFICATION "${PACKAGE_CLASSIFICATION}")
  ENDIF()

  MULTILINE_SET(DOCSTR
    "Enable the package ${PACKAGE_NAME}.  Set to 'ON', 'OFF', or leave"
    " empty to allow for other logic to decide."
    )
  SET_CACHE_ON_OFF_EMPTY( ${PROJECT_NAME}_ENABLE_${PACKAGE_NAME}
    "${PACKAGE_ENABLE}" ${DOCSTR} )

ENDMACRO()


#
# Function that determines if it is okay to allow an implicit package enable
# based on its classification.
#

FUNCTION(PACKAGE_ARCH_IMPLICIT_PACKAGE_ENABLE_IS_ALLOWED PACKAGE_NAME_IN
  IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT
  )
  IF (${PACKAGE_NAME_IN}_CLASSIFICATION STREQUAL PS)
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED TRUE)
  ELSEIF (${PACKAGE_NAME_IN}_CLASSIFICATION STREQUAL SS
    AND ${PROJECT_NAME}_ENABLE_SECONDARY_STABLE_CODE
    )
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED TRUE)
  ELSE()
    SET(IMPLICIT_PACKAGE_ENABLE_ALLOWED FALSE)
  ENDIF()
  SET(${IMPLICIT_PACKAGE_ENABLE_ALLOWED_OUT} ${IMPLICIT_PACKAGE_ENABLE_ALLOWED}
    PARENT_SCOPE )
ENDFUNCTION()


#
# Macro that processes ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS into
# ${PROJECT_NAME}_PACKAGES, ${PROJECT_NAME}_PACKAGE_DIRS, ${PROJECT_NAME}_NUM_PACKAGES,
# ${PROJECT_NAME}_LAST_PACKAGE_IDX, and ${PROJECT_NAME}_REVERSE_PACKAGES.
#
# This macro also sets up the standard package options along with
# default enables/disables.
#

MACRO(PACKAGE_ARCH_PROCESS_PACKAGES_AND_DIRS_LISTS)

  ADVANCED_OPTION(${PROJECT_NAME}_REMOVE_DEFAULT_PACKAGE_DISABLES
    "Removes all default disables from the packages list.  Used for testing etc."
    OFF )

  #
  # Separate out separate lists of package names and directoires
  #

  # Get the total number of packages defined  

  ASSERT_DEFINED(${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  #PRINT_VAR(${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS)
  LIST(LENGTH ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
    ${PROJECT_NAME}_NUM_PACKAGES_AND_FIELDS )
  MATH(EXPR ${PROJECT_NAME}_NUM_PACKAGES
    "${${PROJECT_NAME}_NUM_PACKAGES_AND_FIELDS}/${PLH_NUM_FIELDS_PER_PACKAGE}")
  #PRINT_VAR(${PROJECT_NAME}_NUM_PACKAGES)
  MATH(EXPR ${PROJECT_NAME}_LAST_PACKAGE_IDX "${${PROJECT_NAME}_NUM_PACKAGES}-1")
  
  # Process each of the packages defined

  #PRINT_VAR(APPEND_TO_PACKAGES_LIST)
  IF (NOT APPEND_TO_PACKAGES_LIST)
    #MESSAGE("Wiping clean list of packages")
    SET(${PROJECT_NAME}_PACKAGES)
    SET(${PROJECT_NAME}_PACKAGE_DIRS)
  ENDIF()

  FOREACH(PACKAGE_IDX RANGE ${${PROJECT_NAME}_LAST_PACKAGE_IDX})

    #PRINT_VAR(${PROJECT_NAME}_PACKAGES)

    MATH(EXPR PACKAGE_NAME_IDX "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+0")
    LIST(GET ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_NAME_IDX} PACKAGE )

    MATH(EXPR PACKAGE_DIR_IDX
      "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_DIR_OFFSET}")
    LIST(GET ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_DIR_IDX} PACKAGE_DIR )

    MATH(EXPR PACKAGE_CLASSIFICATION_IDX
      "${PACKAGE_IDX}*${PLH_NUM_FIELDS_PER_PACKAGE}+${PLH_NUM_PACKAGE_CLASSIFICATION_OFFSET}")
    LIST(GET ${PROJECT_NAME}_PACKAGES_AND_DIRS_AND_CLASSIFICATIONS
      ${PACKAGE_CLASSIFICATION_IDX} PACKAGE_CLASSIFICATION )

    SET(PACKAGE_ABS_DIR "${${PROJECT_NAME}_DEPS_HOME_DIR}/packages/${PACKAGE_DIR}")

    IF (EXISTS ${PACKAGE_ABS_DIR})
      SET(PACKAGE_EXISTS TRUE)
    ELSE()
      SET(PACKAGE_EXISTS FALSE)
    ENDIF()

    IF (${PROJECT_NAME}_ASSERT_MISSING_PACKAGES AND NOT PACKAGE_EXISTS)
      MESSAGE(
        "\n***"
        "\n*** Error, the package ${PACKAGE} directory ${PACKAGE_ABS_DIR} does not exist!"
        "\n***\n" )
      MESSAGE(FATAL_ERROR "Stopping due to above error!")
    ENDIF()

    IF (PACKAGE_EXISTS OR ${PROJECT_NAME}_IGNORE_PACKAGE_EXISTS_CHECK)
      LIST(APPEND ${PROJECT_NAME}_PACKAGES ${PACKAGE})
      LIST(APPEND ${PROJECT_NAME}_PACKAGE_DIRS ${PACKAGE_DIR})
      PACKAGE_ARCH_INSERT_STANDARD_PACKAGE_OPTIONS(${PACKAGE} ${PACKAGE_CLASSIFICATION})
      SET(${PACKAGE}_PARENT_PACKAGE "")
    ELSE()
      IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
        MESSAGE(
          "\n***"
          "\n*** WARNING: Excluding package ${PACKAGE} because ${PACKAGE_ABS_DIR}"
            " does not exist!"
          "\n***\n" )
      ENDIF()
    ENDIF()

    #PRINT_VAR(${PROJECT_NAME}_PACKAGES)

  ENDFOREACH()

  # Get the actual number of packages that actually exist

  LIST(LENGTH ${PROJECT_NAME}_PACKAGES ${PROJECT_NAME}_NUM_PACKAGES )
  PRINT_VAR(${PROJECT_NAME}_NUM_PACKAGES)
  MATH(EXPR ${PROJECT_NAME}_LAST_PACKAGE_IDX "${${PROJECT_NAME}_NUM_PACKAGES}-1")

  # Print the final set of packages in debug mode
  
  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    PRINT_VAR(${PROJECT_NAME}_PACKAGES)
    PRINT_VAR(${PROJECT_NAME}_PACKAGE_DIRS)
  ENDIF()
  
  # Create a reverse list for later use
  
  SET(${PROJECT_NAME}_REVERSE_PACKAGES ${${PROJECT_NAME}_PACKAGES})
  LIST(REVERSE ${PROJECT_NAME}_REVERSE_PACKAGES)

ENDMACRO()
