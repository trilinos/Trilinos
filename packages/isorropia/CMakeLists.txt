

#
# A) Define the package
#

TRIBITS_PACKAGE(Isorropia)

#
# B) Set up package-specific options
#

# ToDo: Add these as you need them!
TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

ADD_SUBDIRECTORY(utils)
# ToDo: Only build the above library if we are building tests/examples!

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
