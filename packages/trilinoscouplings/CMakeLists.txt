

#
# A) Define the package
#

TRIBITS_PACKAGE(TrilinosCouplings  DISABLE_CIRCULAR_REF_DETECTION_FAILURE)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Newp_swahili
  HAVE_NEWP_SWAHILI
  "Enable Swahili support."
  OFF )

  SET(TrilinosCouplings_HAVE_ML_NOX_DEFAULT OFF)

TRIBITS_ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ML_Nox
  HAVE_ML_NOX
  "Enable ML with NOX support."
  ${TrilinosCouplings_HAVE_ML_NOX_DEFAULT} )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
