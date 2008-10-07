
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

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

TRILINOS_PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
