INCLUDE(PackageMacros)

#
# A) Define the package
#

PACKAGE(Amesos2)

SET(${PACKAGE_NAME}_VERSION "0.1d")

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_KLU2
  HAVE_AMESOS2_KLU2
  "Enable KLU2 functionality."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLU
  HAVE_AMESOS2_SUPERLU
  "Enable Superlu functionality."
  ON )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SuperLUMT
  HAVE_AMESOS2_SUPERLUMT
  "Enable SuperLU_MT functionality."
  OFF )

PACKAGE_ADD_DEBUG_OPTION()

PACKAGE_ADD_EXPLICIT_INSTANTIATION_OPTION()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
