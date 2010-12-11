
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  pseudo_cg_hb_indefinite
  SOURCES test_pseudo_cg_indefinite.cpp   
  COMM serial mpi
  ARGS
    "--debug --verbose --tol=1e-4"
  STANDARD_PASS_OUTPUT
  )

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
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/BlockCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS bl_cg_hb
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF(${PACKAGE_NAME}_ENABLE_Ifpack)

    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      bl_pcg_hb
      SOURCES test_bl_pcg_hb.cpp createEpetraProblem.cpp
      COMM serial mpi
      ARGS
        "--verbose --filename=bcsstk14.hb"
      STANDARD_PASS_OUTPUT 
      )

  ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
