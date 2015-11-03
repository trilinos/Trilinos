

#
# A) Define the package
#

TRIBITS_PACKAGE(Teko)

#
# B) Set up package-specific options
#

# mfh 02 Nov 2015: Teko currently wants GlobalOrdinal (GO) = int by
# default, and only wants to use GO = long long if
# Teko_USE_LONGLONG_GO (see below) is set.  I want to preserve this
# behavior, but allow Tpetra to disable GO = int.

IF (DEFINED Tpetra_INST_INT_INT)
  IF (Tpetra_INST_INT_INT)
    SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT ON)
  ELSE ()
    SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT OFF)
  ENDIF ()
ELSE () # NOT DEFINED Tpetra_INST_INT_INT
  IF (DEFINED Tpetra_INST_INT_LONG_LONG)
    IF (Tpetra_INST_INT_LONG_LONG)
      SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT ON)
    ELSE ()
      # FIXME (mfh 02 Nov 2015) If Tpetra enables neither
      # GlobalOrdinal = long long nor GlobalOrdinal = int,
      # I can't preserve Teko's current behavior.  The right
      # thing for Teko to do would be to have logic to use
      # any enabled Tpetra GlobalOrdinal type.
      SET(${PACKAGE_NAME}_USE_LONGLONG_GO_DEFAULT OFF)
    ENDIF ()
  ENDIF ()
ENDIF () # DEFINED Tpetra_INST_INT_INT

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
