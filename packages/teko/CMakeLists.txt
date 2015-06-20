

#
# A) Define the package
#

TRIBITS_PACKAGE(Teko)

#
# B) Set up package-specific options
#

#TRIBITS_ADD_OPTION_AND_DEFINE(Teko_ENABLE_BLOCK_SMOOTHERS
#  Teko_ENABLE_ML_SMOOTHERS
#  "Enable ml block smoothers."
#  OFF )

# Enable ML block Smoothers
IF(ML_ENABLE_TekoSmoothers)
   SET(Teko_ENABLE_ML_SMOOTHERS ON)
ELSE()
   SET(Teko_ENABLE_ML_SMOOTHERS OFF)
ENDIF()

# enable development mode
IF(${PROJECT_NAME}_ENABLE_DEVELOPMENT_MODE)
   SET(Teko_ENABLE_DEV_MODE ON)
ELSE()
   SET(Teko_ENABLE_DEV_MODE OFF)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(tests)
TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)

TRIBITS_EXCLUDE_FILES(
  tests/data/nsjac.mm
  tests/data/nsjac_1.mm
  src/mlutils
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
