
INCLUDE(PackageMacros)
INCLUDE(AddOptionAndDefine)

#
# A) Define the package
#

PACKAGE(TrilinosCouplings)

#
# B) Set up package-specific options
#

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_Newp_swahili
  HAVE_NEWP_SWAHILI
  "Enable Swahili support."
  OFF )

IF(Trilinos_ENABLE_ML AND Trilinos_ENABLE_NOX)
  SET(TrilinosCouplings_HAVE_ML_NOX_DEFAULT ON)
ELSE()
  SET(TrilinosCouplings_HAVE_ML_NOX_DEFAULT OFF)
ENDIF()

ADD_OPTION_AND_DEFINE(${PACKAGE_NAME}_ENABLE_ML_Nox
  HAVE_ML_NOX
  "Enable ML with NOX support."
  ${TrilinosCouplings_HAVE_ML_NOX_DEFAULT} )

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

#PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(examples)

#
# D) Do standard postprocessing
#

PACKAGE_POSTPROCESS()
