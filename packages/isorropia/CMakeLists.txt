
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Isorropia)

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

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
