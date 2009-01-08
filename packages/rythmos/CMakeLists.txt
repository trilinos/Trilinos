
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Rythmos)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
