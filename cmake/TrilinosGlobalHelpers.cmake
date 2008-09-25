
#
# Macro that sets up standard user options for each package
#

MACRO(TRILINOS_INSERT_STANDARD_PACAKGE_VARIABLES PACAKGE)

  #
  # Global user-level variables that define what gets enabled
  #

  SET( Trilinos_ENABLE_${PACKAGE} "" CACHE BOOL
    "Enable the ${PACKAGE} package.")
  SET( ${PACKAGE}_ENABLE_TESTS "" CACHE BOOL
    "Build ${PACKAGE} tests." )
  SET( ${PACKAGE}_ENABLE_EXAMPLES "" CACHE BOOL
    "Build ${PACKAGE} examples." )

  #
  # Global internal variables that link everything together
  #

  SET(${PACKAGE}_DEP_PACKAGES "" CACHE INTERNAL "${PACKAGE} dependent packages")
  SET(${PACKAGE}_INCLUDE_DIRS "" CACHE INTERNAL "${PACKAGE} include directories")
  SET(${PACKAGE}_LIBRARY_DIRS "" CACHE INTERNAL "${PACKAGE} library directories")
  SET(${PACKAGE}_LIBRARIES "" CACHE INTERNAL "${PACKAGE} libraries")

ENDMACRO()


#
# Set an individual pacakge variable based on the global value
#

MACRO(TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLE TRILINOS_VAR PACKAGE_VAR)

  IF (VERBOSE_CONFIGURE)
    MESSAGE("")
    MESSAGE("TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLE:")
    MESSAGE(STATUS "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
    MESSAGE(STATUS "${TRILINOS_VAR} = ${${TRILINOS_VAR}}")
  ENDIF()

  IF (${PACKAGE_VAR} STREQUAL "")
    IF (${TRILINOS_VAR} STREQUAL "ON")
      #MESSAGE(STATUS "ON")
      SET(${PACKAGE_VAR} ON)
    ELSEIF (TRILINOS_VAR STREQUAL "OFF")
      #MESSAGE(STATUS "OFF")
      SET(${PACKAGE_VAR} OFF)
    ELSE()
      #MESSAGE(STATUS "ELSE")
      # Otherwise, we will leave it up the the individual package
      # to decide?
    ENDIF()
  ELSE()
    #MESSAGE(STATUS "PACAKGE_VAR NOT DEFAULT")
  ENDIF()

  IF (VERBOSE_CONFIGURE)
    MESSAGE(STATUS "${PACKAGE_VAR} = ${${PACKAGE_VAR}}")
  ENDIF()

ENDMACRO()


#
# Macro that enables and disables individual package options based on global
# options.
#
# Depends on the global variables:
#
#   Trilinos_ENABLE_ALL_PACKAGES
#   Trilinos_ENABLE_TESTS
#   Trilinos_ENABLE_EXAMPLES
#
# Which are set by the user
#

MACRO(TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLES PACAKGE)

  TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLE(
    Trilinos_ENABLE_ALL_PACKAGES Trilinos_ENABLE_${PACKAGE} )

  TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLE(
    Trilinos_ENABLE_TESTS ${PACKAGE}_ENABLE_TESTS )

  TRILINOS_POSTPROCESS_STANDARD_PACAKGE_VARIABLE(
    Trilinos_ENABLE_EXAMPLES ${PACKAGE}_ENABLE_EXAMPLES )

ENDMACRO()
