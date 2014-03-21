#
# A) Define the package
#
TRIBITS_PACKAGE( SimpleCxx  ENABLE_SHADOWING_WARNINGS  CLEANED )

#
# B) Platform-specific checks
#
INCLUDE(CheckFor__int64)
CHECK_FOR___INT64(HAVE_SIMPLECXX___INT64)

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
