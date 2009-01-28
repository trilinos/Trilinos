
INCLUDE(PackageAddExecutable)
INCLUDE(PackageAddTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    bl_cg_hb
    SOURCES test_bl_cg_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=bcsstk14.hb"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    resolve_cg_hb
    SOURCES test_resolve_cg_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    ARGS
      "--verbose --filename=bcsstk14.hb"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/example/BlockCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS bl_cg_hb
    )

  ASSERT_DEFINED(Trilinos_ENABLE_Ifpack)
  IF(Trilinos_ENABLE_Ifpack)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_pcg_hb
      SOURCES test_bl_pcg_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=bcsstk14.hb"
      STANDARD_PASS_OUTPUT 
      )

  ENDIF(Trilinos_ENABLE_Ifpack)

  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_cg_complex_hb
      SOURCES test_bl_cg_complex_hb.cpp createEpetraProblem.cpp
      ARGS
        "--verbose --filename=mhd1280b.cua"
      )

    PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestBlockCGComplexFiles
      SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/test/MVOPTester
      SOURCE_FILES mhd1280b.cua
      EXEDEPS bl_cg_hb
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
