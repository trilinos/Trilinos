
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(Pliris)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ZCPLX
  ZCPLX
  "Enable double precision complex functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SCPLX
  SCPLX
  "Enable single precision complex functionality."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DREAL
  DREAL
  "Set reals to double precision."
  OFF )

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SREAL
  SREAL
  "Set reals to single precision."
  OFF )

#
# C) Add the libraries, tests, and examples
#

#If no options are specified set DREAL as the default
IF(NOT ${PACKAGE_NAME}_ENABLE_ZCPLX AND
   NOT ${PACKAGE_NAME}_ENABLE_SCPLX AND
   NOT ${PACKAGE_NAME}_ENABLE_SREAL AND
   NOT ${PACKAGE_NAME}_ENABLE_DREAL)
   
  SET(${PACKAGE_NAME}_ENABLE_DREAL ON)
ENDIF()

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
