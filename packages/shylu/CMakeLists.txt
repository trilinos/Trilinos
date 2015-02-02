

#
# A) Define the package
#

TRIBITS_PACKAGE_DECL(ShyLU)

#
# B) Set up package-specific options
#

#TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
#  HAVE_METIS
#  "Metis is present"
#  OFF )


#
# C) Add the libraries, tests, and examples
#

## ADD_SUBDIRECTORY(src)

## TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

TRIBITS_PROCESS_SUBPACKAGES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
