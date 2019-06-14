

#
# A) Define the package
#

TRIBITS_PACKAGE(Teko)

#
# B) Set up package-specific options
#

# FIXME (mfh 02 Nov 2015, 04 Nov 2015) The right thing for Teko to do
# would be to have logic to use any enabled Tpetra GlobalOrdinal type.
# For now, we just error out if neither long long nor int are enabled.

IF (DEFINED Tpetra_INST_INT_LONG_LONG)
  IF (Tpetra_INST_INT_LONG_LONG)
    SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT ON)
  ELSE ()
    IF (DEFINED Tpetra_INST_INT_INT)
      IF (Tpetra_INST_INT_INT)
	SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT OFF)
      ELSE ()
	MESSAGE (FATAL_ERROR "Teko currently requires that Tpetra enable either GlobalOrdinal (GO) = long long or GO = int.  Please set Tpetra_INST_INT_LONG_LONG=ON to enable GO = long long, or Tpetra_INST_INT_INT=ON to enable GO = int.  If you would like to use other GO types, please contact the Teko developers.")
      ENDIF () # Tpetra_INST_INT_INT
    ELSE () # NOT DEFINED Tpetra_INST_INT_INT
      MESSAGE (FATAL_ERROR "Teko currently requires that Tpetra enable either GlobalOrdinal (GO) = long long or GO = int.  Please set Tpetra_INST_INT_LONG_LONG=ON to enable GO = long long, or Tpetra_INST_INT_INT=ON to enable GO = int.  If you would like to use other GO types, please contact the Teko developers.")
    ENDIF ()
  ENDIF () # Tpetra_INST_INT_LONG_LONG
ELSE () # NOT DEFINED Tpetra_INST_INT_LONG_LONG
  IF (NOT DEFINED Tpetra_INST_INT_INT)
    MESSAGE (FATAL_ERROR "Teko currently requires that Tpetra enable either GlobalOrdinal (GO) = long long or GO = int.  Please set Tpetra_INST_INT_LONG_LONG=ON to enable GO = long long, or Tpetra_INST_INT_INT=ON to enable GO = int.  If you would like to use other GO types, please contact the Teko developers.")
  ELSEIF (NOT Tpetra_INST_INT_INT)
    MESSAGE (FATAL_ERROR "Teko currently requires that Tpetra enable either GlobalOrdinal (GO) = long long or GO = int.  Please set Tpetra_INST_INT_LONG_LONG=ON to enable GO = long long, or Tpetra_INST_INT_INT=ON to enable GO = int.  If you would like to use other GO types, please contact the Teko developers.")
  ENDIF ()
  SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT OFF)
ENDIF ()

ASSERT_DEFINED(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT)
TRIBITS_ADD_OPTION_AND_DEFINE(
  ${PACKAGE_NAME}_USE_LONGLONG_GO
  HAVE_${PACKAGE_NAME}_USE_LONGLONG_GO
  "Use long long for global ordinal"
  ${${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT}
  )

ASSERT_DEFINED(${PACKAGE_NAME}_USE_LONGLONG_GO)
IF (${PACKAGE_NAME}_USE_LONGLONG_GO)
  IF (DEFINED Tpetra_INST_INT_LONG_LONG)
    IF (NOT Tpetra_INST_INT_LONG_LONG)
      MESSAGE(FATAL_ERROR "Teko: If GO=long long is disabled in Tpetra, you may not set ${PACKAGE_NAME}_USE_LONGLONG_GO = ON.  First, set Tpetra_INST_INT_LONG_LONG=ON to enable GO=long long in Tpetra.")
    ENDIF ()
  ELSE () # NOT DEFINED Tpetra_INST_INT_LONG_LONG
    MESSAGE(FATAL_ERROR "Teko: If GO=long long is disabled in Tpetra, you may not set ${PACKAGE_NAME}_USE_LONGLONG_GO = ON.  First, set Tpetra_INST_INT_LONG_LONG=ON to enable GO=long long in Tpetra.")      
  ENDIF ()
ENDIF ()


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
