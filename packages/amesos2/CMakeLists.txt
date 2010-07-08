
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Amesos2)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_EXPERIMENTAL
  HAVE_EXPERIMENTAL
  "Enable experimental features"
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_BUILD_BTF
  HAVE_BTF
  "Enable BTF functionality"
  OFF )

ADD_OPTION_AND_DEFINE(Trilinos_ENABLE_PyTrilinos
  HAVE_PYTRILINOS
  "Enable certain python-compatible code only if PyTrilinos is enabled."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Superlu
  HAVE_AMESOS_SUPERLU
  "Enable Superlu functionality."
  ON )

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
