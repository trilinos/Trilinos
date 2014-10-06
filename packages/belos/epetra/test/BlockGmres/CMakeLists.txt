

ASSERT_DEFINED(Belos_ENABLE_Triutils)
IF (Belos_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    bl_gmres_diag
    SOURCES test_bl_gmres_diag.cpp createEpetraProblem.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    bl_gmres_hb
    SOURCES test_bl_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    resolve_gmres_hb
    SOURCES test_resolve_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_gmres_hb
    SOURCES test_pseudo_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS "--verbose --filename=orsirr1.hb"
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/epetra/example/BlockGmres
    SOURCE_FILES orsirr1.hb fidap036.hb
    EXEDEPS bl_gmres_diag
    )

  ASSERT_DEFINED(Belos_ENABLE_Ifpack)
  IF(Belos_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bl_pgmres_hb
      SOURCES test_bl_pgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=orsirr1.hb"
        "--verbose --filename=orsirr1.hb --num-rhs=1 --blocksize=2"
        "--verbose --filename=orsirr1.hb --pseudo"
      STANDARD_PASS_OUTPUT
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      pseudo_pgmres_hb
      SOURCES test_pseudo_pgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bl_fgmres_hb
      SOURCES test_bl_fgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT
      )

  ENDIF(Belos_ENABLE_Ifpack)

ENDIF(Belos_ENABLE_Triutils)
