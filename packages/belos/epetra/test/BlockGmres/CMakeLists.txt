
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  bl_gmres_diag
  SOURCES test_bl_gmres_diag.cpp 
  COMM serial mpi
  STANDARD_PASS_OUTPUT
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    bl_gmres_hb
    SOURCES test_bl_gmres_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    resolve_gmres_hb
    SOURCES test_resolve_gmres_hb.cpp 
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    pseudo_gmres_hb
    SOURCES test_pseudo_gmres_hb.cpp 
    COMM serial mpi
    ARGS "--verbose --filename=orsirr1.hb"
    STANDARD_PASS_OUTPUT 
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES orsirr1.hb 
    EXEDEPS bl_gmres_hb
  )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bl_pgmres_hb
      SOURCES test_bl_pgmres_hb.cpp 
      COMM serial mpi
      ARGS
        "--verbose --filename=orsirr1.hb"
        "--verbose --filename=orsirr1.hb --num-rhs=1 --blocksize=2"
        "--verbose --filename=orsirr1.hb --pseudo"
        "--verbose --filename=nos1.rsa --num-rhs=8 --blocksize=2 --ortho=ICGS"
        "--verbose --filename=nos1.rsa --num-rhs=8 --blocksize=2 --ortho=IMGS"
        "--verbose --filename=nos1.rsa --num-rhs=8 --blocksize=2 --ortho=DGKS"
      STANDARD_PASS_OUTPUT 
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      pseudo_pgmres_hb
      SOURCES test_pseudo_pgmres_hb.cpp 
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT 
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      bl_fgmres_hb
      SOURCES test_bl_fgmres_hb.cpp 
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT 
      )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyTestBlockPGmresFiles
      SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
      SOURCE_FILES nos1.rsa
      EXEDEPS bl_pgmres_hb
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
