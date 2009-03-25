
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(AztecOO)

#
# B) Set up package-specific options
#
ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_AZLU
  HAVE_AZLU
  "Enables the azlu preconditioner.  Requires the y12m TPL."
  "${${PACKAGE_NAME}_ENABLE_y12m}")

IF (NOT ${PACKAGE_NAME}_ENABLE_y12m AND ${PACKAGE_NAME}_ENABLE_AZLU)
  MESSAGE(FATAL_ERROR "You have enabled the AZLU support in AztecOO with the flag \"AztecOO_ENABLE_AZLU\", but this requires the y12m TPL which is currently disabled.  Please enable this TPL or disable AZLU support.")
ENDIF()

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
