#
# A) Define the package
#

TRIBITS_PACKAGE(  MixedLanguage  ENABLE_SHADOWING_WARNINGS  CLEANED )

#
# B) Platform-specific checks
#

# ToDo: Give an example

#
# C) Set up package-specific options
#

TRIBITS_ADD_DEBUG_OPTION()

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# D) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# E) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
