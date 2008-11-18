
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Tpetra)

#
# B) Set up package-specific options
#

TRILINOS_ADD_OPTION(${PACKAGE_NAME}_ENABLE_TBB
  HAVE_TPETRA_TBB
  "Enable Tpetra's use of Intel's Thread Building Blocks library"
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

#TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
