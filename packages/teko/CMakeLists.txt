

#
# A) Define the package
#

TRIBITS_PACKAGE(Teko)

#
# B) Set up package-specific options
#

TRIBITS_ADD_OPTION_AND_DEFINE(
  Ifpack2_Trilinos
  HAVE_${PACKAGE_NAME}_ENABLE_Ifpack2
  "Enable Ifpack2 preconditioners"
  OFF
  )

#TRIBITS_ADD_OPTION_AND_DEFINE(Teko_ENABLE_BLOCK_SMOOTHERS
#  Teko_ENABLE_ML_SMOOTHERS
#  "Enable ml block smoothers."
#  OFF )

IF(${PROJECT_NAME}_ENABLE_Epetra AND DEFINED Teko_ENABLE_Epetra)
   GLOBAL_SET(TEKO_HAVE_EPETRA ${Teko_ENABLE_Epetra})
ELSE()
   GLOBAL_SET(TEKO_HAVE_EPETRA ${${PROJECT_NAME}_ENABLE_Epetra})
ENDIF()

MESSAGE(STATUS "TEKO_HAVE_EPETRA = ${TEKO_HAVE_EPETRA}")

IF(${PACKAGE_NAME}_ENABLE_MueLu)
   GLOBAL_SET(HAVE_MUELU ON)
ELSE()
   GLOBAL_SET(HAVE_MUELU OFF)
ENDIF()

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

# enable development mode
IF(${PROJECT_NAME}_ENABLE_DEBUG_MSGS)
   SET(Teko_DEBUG_OFF OFF)
ELSE()
   SET(Teko_DEBUG_OFF ON)
ENDIF()

#
# C) Add the libraries, tests, and examples
#

ADD_SUBDIRECTORY(src)

TRIBITS_ADD_TEST_DIRECTORIES(tests)
IF(TEKO_HAVE_EPETRA)
   TRIBITS_ADD_EXAMPLE_DIRECTORIES(examples)
ENDIF()

TRIBITS_EXCLUDE_FILES(
  tests/data/nsjac.mm
  tests/data/nsjac_1.mm
  src/mlutils
  )

#
# D) Do standard postprocessing
#

TRIBITS_PACKAGE_POSTPROCESS()
