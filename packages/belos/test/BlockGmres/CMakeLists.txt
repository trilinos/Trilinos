
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    bl_gmres_diag
    SOURCES test_bl_gmres_diag.cpp createEpetraProblem.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    bl_gmres_hb
    SOURCES test_bl_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    resolve_gmres_hb
    SOURCES test_resolve_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=orsirr1.hb"
      "--verbose --filename=orsirr1.hb --pseudo"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    pseudo_gmres_hb
    SOURCES test_pseudo_gmres_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS "--verbose --filename=orsirr1.hb"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/example/BlockGmres
    SOURCE_FILES orsirr1.hb fidap036.hb
    EXEDEPS bl_gmres_diag
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_pgmres_hb
      SOURCES test_bl_pgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=orsirr1.hb"
        "--verbose --filename=orsirr1.hb --pseudo"
      STANDARD_PASS_OUTPUT 
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      pseudo_pgmres_hb
      SOURCES test_pseudo_pgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT 
      )

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_fgmres_hb
      SOURCES test_bl_fgmres_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS "--verbose --filename=orsirr1.hb"
      STANDARD_PASS_OUTPUT 
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_gmres_complex_hb
      SOURCES test_bl_gmres_complex_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=mhd1280b.cua"
        "--verbose --filename=mhd1280b.cua --pseudo"
        "--verbose --filename=mhd1280a.cua --subspace-length=200 --tol=1e-4"
      STANDARD_PASS_OUTPUT 
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockGmresComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua mhd1280a.cua
      EXEDEPS bl_gmres_diag
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
