
INCLUDE(TribitsPackageMacros)
INCLUDE(TribitsAddOptionAndDefine)

#
# A) Define the package
#

TRIBITS_PACKAGE(ROL)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DEBUG
  HAVE_ROL_DEBUG
  "Enable a host of runtime debug checking."
  ${Trilinos_ENABLE_DEBUG}
  )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test util)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_AUTOTOOLS_FILES()

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
