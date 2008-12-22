INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Phalanx)

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_EXPLICIT_TEMPLATE_INSTANTIATION
  PHX_ETI
  "Build with explicit template instation enabled."
  OFF )

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_DEBUG
  PHX_DEBUG
  "Build with explicit template instation enabled."
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
