
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Isorropia)

#
# B) Set up package-specific options
#

# ToDo: Add these as you need them!

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(utils)
# ToDo: Only build the above library if we are building tests/examples!

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

TRILINOS_PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
