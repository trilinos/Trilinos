ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Tpetra_BiCGStab_diag
  SOURCES test_bicgstab_diag.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_BiCGStab_hb_test
    SOURCES test_bicgstab_hb.cpp
    COMM serial mpi
    ARGS "--verbose --max-iters=10"
         "--verbose --max-iters=10 --tol=1e-10"
         "--verbose --max-iters=10 --num-rhs=2"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Tpetra_KokkosDense_BiCGStab_hb_test
    SOURCES test_kokkos_bicgstab_hb.cpp
    COMM serial mpi
    ARGS "--verbose --max-iters=10"
         "--verbose --max-iters=10 --tol=1e-10"
         "--verbose --max-iters=10 --num-rhs=2"
    STANDARD_PASS_OUTPUT
    )

IF (Tpetra_INST_COMPLEX_DOUBLE)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Tpetra_bicgstab_complex_hb
      SOURCES test_bicgstab_complex_hb.cpp
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    ASSERT_DEFINED(Anasazi_SOURCE_DIR)
    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTpetraBiCGStabComplexFiles
      SOURCE_DIR ${Anasazi_SOURCE_DIR}/testmatrices
      SOURCE_FILES mhd1280b.cua
      EXEDEPS Tpetra_bicgstab_complex_hb
    )

ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_CopyTestBiCGSTABFiles
    SOURCE_DIR ${${PACKAGE_NAME}_SOURCE_DIR}/testmatrices
    SOURCE_FILES cage4.hb
    EXEDEPS Tpetra_BiCGStab_hb_test Tpetra_KokkosDense_BiCGStab_hb_test
    )

ENDIF()
