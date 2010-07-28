
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Teko)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Isorropia
   Teko_ENABLE_Isorropia
   "Turn on Isorropia related options including probing"
   OFF)

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(tests)
PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
