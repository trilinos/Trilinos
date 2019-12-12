

#
# A) Define the package
#

TRIBITS_PACKAGE(Claps)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SClaps
  HAVE_SCLAPS
  "Enable Sandia Claps interface."
  OFF )

# Claps package is deprecated April 2019; to be removed May 2019
MESSAGE(WARNING "***** WARNING:  Package Claps is deprecated in April 2019; it will be removed in May 2019 *****")

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
