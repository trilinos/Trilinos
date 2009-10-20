
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Teko)

#
# B) Set up package-specific options
#

# Nothing to do here

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(tests)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
