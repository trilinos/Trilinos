
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Teko)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_DEVELOPMENT_MODE
  Teko_ENABLE_DEV_MODE
  "Enable development functionality."
  ON )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(tests)
PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

PACKAGE_ARCH_EXCLUDE_FILES(
  tests/data/nsjac.mm
  tests/data/nsjac_1.mm
  src/mlutils
  )

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
