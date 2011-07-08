INCLUDE(PackageMacros)


#
# Define a subpackage.
#
# Once called, the following varibles are in scope:
#
#   PARENT_PACKAGE_NAME: The name of the parent package
#
#   SUBPACKAGE_NAME: The local name of the subpackage (does not contain the
#   parent package name).
#
#   SUBPACKAGE_FULLNAME: The full project-level name of the subpackage (which
#   includes the parent package name at the beginning).
#
#   PACKAGE_NAME: Inside the subpackage, the same as SUBPACKAGE_FULLNAME.
#
MACRO(SUBPACKAGE SUBPACKAGE_NAME_IN)

  IF (${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE("\nSUBPACKAGE: ${SUBPACKAGE_NAME_IN}")
  ENDIF()

  IF (NOT ${SUBPACKAGE_NAME_IN} STREQUAL ${SUBPACKAGE_NAME})
    MESSAGE(FATAL_ERROR "Error, the package-defined subpackage name"
      " '${SUBPACKAGE_NAME_IN}' is not the same as the subpackage name"
      " '${SUBPACKAGE_NAME}' defined in the parent packages's"
      " Dependencies.cmake file")
  ENDIF()

  # To provide context for various macros
  SET(PACKAGE_NAME ${SUBPACKAGE_FULLNAME})

  SET(PARENT_PACKAGE_SOURCE_DIR "${PACKAGE_SOURCE_DIR}")
  SET(PARENT_PACKAGE_BINARY_DIR "${PACKAGE_BINARY_DIR}")

  # Now override the package-like varibles
  PACKAGE_SET_COMMON_VARS(${SUBPACKAGE_FULLNAME})
  PACKAGE_DEFINE_LINKAGE_VARS(${SUBPACKAGE_FULLNAME})

ENDMACRO()


#
# Postprocess after defining a subpackage
#

MACRO(SUBPACKAGE_POSTPROCESS)
  PACKAGE_POSTPROCESS()
ENDMACRO()
