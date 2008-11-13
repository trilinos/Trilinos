
INCLUDE(TrilinosPackageHelpers)
INCLUDE(Trilinos_Add_Option)

#
# A) Define the package
#

TRILINOS_PACKAGE(Galeri)

#
# B) Set up package-specific options
#

# ...

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#this dir only had "EXTRA_DIST" files in it
#ADD_SUBDIRECTORY(src-fem)  

TRILINOS_PACKAGE_ADD_TEST_DIRECTORIES(test)

TRILINOS_PACKAGE_ADD_EXAMPLE_DIRECTORIES(example example-fem)

#
# D) Do standard postprocessing
#

TRILINOS_PACKAGE_POSTPROCESS()
