
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Teuchos)

SET(${PACKAGE_NAME}_VERSION "1.5d")

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(
  ${PACKAGE_NAME}_ENABLE_COMPLEX
  HAVE_TEUCHOS_COMPLEX
  "Enable Teuchos features for complex arithmetic."
  ON
  )

TRILINOS_ADD_OPTION(
  ${PACKAGE_NAME}_ENABLE_EXTENDED
  HAVE_TEUCHOS_EXTENDED
  "Enable extended Teuchos features."
  ON
  )

TRILINOS_ADD_OPTION(
  ${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_TEUCHOS_DEBUG
  "Enable a host of runtime debug checking."
  OFF
  )

TRILINOS_ADD_OPTION(
  ${PACKAGE_NAME}_ENABLE_GCC_DEMANGLE
  HAVE_TEUCHOS_DEMANGLE
  "Enable Teuchos calls to GCC demangle library function."
  OFF
  )

# 2008/10/02: rabartl: ToDo: Above: I need to add a macro like
# TRILINOS_DECLARE_OPTION(OPTION_NAME DEFINE_NAME DOCUMENTATION)

#
# C) Add the libraries and examples
#

ADD_SUBDIRECTORY(src)

IF(${PACKAGE_NAME}_ENABLE_TESTS)
  ADD_SUBDIRECTORY(test)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_EXAMPLES)
  ADD_SUBDIRECTORY(example)
ENDIF()

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
