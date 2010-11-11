
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(KrylovSupport)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SOMETHING
#  HAVE_${PACKAGE_NAME}_SOMETHING
#  "Enable some functionality."
#  ON )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
