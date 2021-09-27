

#
# A) Forward delcare the package so that certain options are also defined for
# subpackages
#

TRIBITS_PACKAGE_DECL(Thyra ENABLE_SHADOWING_WARNINGS)

#
# B) Define the common options for Thyra first so they can be used by
# subpackages as well.
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_ENABLE_TEUCHOS_TIME_MONITOR_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_EXPLICIT_INSTANTIATION_OPTION()

IF (Trilinos_ENABLE_FLOAT  AND  Teuchos_ENABLE_FLOAT)
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_FLOAT
  HAVE_THYRA_FLOAT
  "Enable Thyra features for float data type."
  ${${PACKAGE_NAME}_ENABLE_FLOAT_DEFAULT} )

IF (Trilinos_ENABLE_COMPLEX_DOUBLE  AND  Teuchos_ENABLE_COMPLEX)
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT  ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT OFF)
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_THYRA_COMPLEX
  "Enable Thyra features for complex data types."
  ${${PACKAGE_NAME}_ENABLE_COMPLEX_DEFAULT} )

# ToDo: Above, split ${PACKAGE_NAME}_ENABLE_COMPLEX into
# ${PACKAGE_NAME}_ENABLE_COMPLEX_FLOAT and
# ${PACKAGE_NAME}_ENABLE_COMPLEX_DOUBLE.  See Trilinos GitHub #362.

IF (MSVC)
  # See Trilinos bug 5785
  SET(${PACKAGE_NAME}_ENABLE_ME_POLYNOMIAL_DEFAULT OFF)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_ME_POLYNOMIAL_DEFAULT
    ${${PROJECT_NAME}_ENABLE_SECONDARY_TESTED_CODE} )
ENDIF()

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_ME_POLYNOMIAL
  HAVE_THYRA_ME_POLYNOMIAL
  "Enable Thyra ModelEvalutor Polynomial InArgs and OutArgs."
  ${${PACKAGE_NAME}_ENABLE_ME_POLYNOMIAL_DEFAULT} )


# Allow for disabling Hessian interface objects in Thyra ModelEvaluator
# and OutArgs. Doing so can significantly reduce overhead (memory footprint
# and computational cost) when large numbers of responses and parameters are
# used and Hessian derivatives are not required.

TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_ENABLE_HESSIAN_SUPPORT
  ${PACKAGE_NAME}_BUILD_HESSIAN_SUPPORT
  "Enable building of Hessian interfaces for Thyra ModelEvaluator."
  ON )

IF(${PACKAGE_NAME}_BUILD_HESSIAN_SUPPORT)
  MESSAGE("-- Thyra ModelEvaluator Hessian interfaces On")
ELSE()
  MESSAGE("-- Thyra ModelEvaluator Hessian interfaces Off")
ENDIF()

#
# C) Process the subpackages for Thyra
#

TRIBITS_PROCESS_SUBPACKAGES()

#
# D) If Thyra itself is enabled, process the Thyra package
#

TRIBITS_PACKAGE_DEF()

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

TRIBITS_PACKAGE_POSTPROCESS()
