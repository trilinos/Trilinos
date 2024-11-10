

#
# A) Define the package
#

TRIBITS_PACKAGE(Pliris)

#
# B) Set up package-specific options
#

TRIBITS_ADD_SHOW_DEPRECATED_WARNINGS_OPTION()

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ZCPLX
  ZCPLX
  "Enable double precision complex functionality."
  OFF NONCACHE)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SCPLX
  SCPLX
  "Enable single precision complex functionality."
  OFF NONCACHE)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_SREAL
  SREAL
  "Set reals to single precision."
  OFF NONCACHE)

#
# C) Add the libraries, tests, and examples
#

#If no options are specified set DREAL as the default
IF(NOT ${PACKAGE_NAME}_ENABLE_ZCPLX AND
   NOT ${PACKAGE_NAME}_ENABLE_SCPLX AND
   NOT ${PACKAGE_NAME}_ENABLE_SREAL
  )
  SET(${PACKAGE_NAME}_ENABLE_DREAL_DEFAULT ON)
ELSE()
  SET(${PACKAGE_NAME}_ENABLE_DREAL_DEFAULT OFF)
ENDIF()
TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_DREAL
  DREAL
  "Set reals to double precision."
  ${${PACKAGE_NAME}_ENABLE_DREAL_DEFAULT} NONCACHE )

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(test)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

#
# Exclude files for source package.
#

TRIBITS_EXCLUDE_FILES(
  doc/matrix_.*.gif
  src/Make..*
  src/clean_code.h
  src/init..*
  src/malloc.c
  src/my_srand48.c
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
