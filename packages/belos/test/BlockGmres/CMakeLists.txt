
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils)
    
  ASSERT_DEFINED(Teuchos_ENABLE_COMPLEX)
  IF(Teuchos_ENABLE_COMPLEX)

    INCLUDE_DIRECTORIES(../MVOPTester)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_gmres_complex_hb
      SOURCES test_bl_gmres_complex_hb.cpp 
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
      EXEDEPS bl_gmres_complex_hb
      )

  ENDIF(Teuchos_ENABLE_COMPLEX)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
