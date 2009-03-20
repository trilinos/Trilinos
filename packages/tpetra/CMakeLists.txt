
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE( Tpetra ENABLE_SHADOWING_WARNINGS )
# CLEANED

# 2009/03/20: rabartl: I removed 'CLEANED' above because there were warnings
# and the code did not build.  Once we have made the switch to CMake as the
# default build system and we have a pre-checkin testing tool we can enable
# these again.

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_TBB
  HAVE_TPETRA_TBB
  "Enable Tpetra's use of Intel's Thread Building Blocks library"
  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_TEST_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
