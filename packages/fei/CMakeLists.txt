
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(FEI)

#
# B) Set up package-specific options
#

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(base)
ADD_SUBDIRECTORY(support-Trilinos)
ADD_SUBDIRECTORY(test_utils)
ADD_SUBDIRECTORY(utest_cases)


PACKAGE_ADD_TEST_DIRECTORIES(test)
PACKAGE_ADD_TEST_DIRECTORIES(utest)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
