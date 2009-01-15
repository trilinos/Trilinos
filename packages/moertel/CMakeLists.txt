
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Moertel)

#
# B) Set up package-specific options
#

#ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Experimental
#  HAVE_KOKKOS_EXPERIMENTAL
#  "Enable experimental."
#  OFF )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
