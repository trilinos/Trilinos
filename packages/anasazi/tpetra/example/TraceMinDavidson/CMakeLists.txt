# FIXME (mfh 16 Jan 2016) The example doesn't currently build
# when TSQR is enabled.  This is because the test uses a custom
# specialization of MultiVecTraits that does not provide the
# 'tsqr_adaptor_type' typedef.  For now, I'm disabling the tests
# when TSQR is enabled.

ASSERT_DEFINED(Anasazi_ENABLE_TSQR)
IF(NOT Anasazi_ENABLE_TSQR)

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_Gen_example
  SOURCES TraceMinDavidsonGeneralizedEx.cpp
  COMM serial mpi
)

# FIXME (mfh 22 Oct 2015) This example requires GlobalOrdinal =
# LocalOrdinal.  It should be easy to fix this, but for now, we
# disable it if GlobalOrdinal = LocalOrdinal = int is disabled.
IF (Tpetra_INST_INT_INT)
  TRIBITS_ADD_EXECUTABLE(
    Tpetra_TD_Fiedler_example
    SOURCES TraceMinDavidsonLaplacianEx.cpp
    COMM serial mpi
    )
ENDIF ()

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_Trans_example
  SOURCES TraceMinDavidsonSpecTransEx.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Tpetra_TD_UserOp_example
  SOURCES TraceMinDavidsonUserOpEx.cpp
  COMM serial mpi
)

ENDIF()
