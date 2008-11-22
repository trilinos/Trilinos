
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Ifpack)

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_METIS
  HAVE_IFPACK_METIS
  "Enable support for Metis in Ifpack."
  OFF )

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
