INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Piro DISABLE_STRONG_WARNINGS)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)
PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
