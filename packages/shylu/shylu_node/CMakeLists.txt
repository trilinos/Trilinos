#
# A) Define the package
#

TRIBITS_PACKAGE_DECL(ShyLU_Node)

#
# C) Add the libraries, tests, and examples
#

## ADD_SUBDIRECTORY(src)

## TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PROCESS_SUBPACKAGES()

TRIBITS_PACKAGE_DEF()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
