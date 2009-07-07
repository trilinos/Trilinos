
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Ifpack)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_METIS
  HAVE_IFPACK_METIS
  "Enable support for Metis in Ifpack."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_USING_HYPRE
  HAVE_HYPRE
  "Enable Ifpack/Hypre adapters"
  OFF )


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
