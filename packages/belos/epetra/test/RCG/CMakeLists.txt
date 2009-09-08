
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
IF (${PACKAGE_NAME}_ENABLE_Triutils AND ${PACKAGE_NAME}_ENABLE_EpetraExt)
    
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    rcg_hb
    SOURCES test_rcg_hb.cpp createEpetraProblem.cpp
    COMM serial mpi
    XHOST gabriel.sandia.gov
    ARGS
      "--verbose --tol=1e-6 --filename=bcsstk14.hb --num-rhs=10 --max-subspace=100 --recycle=10 --max-iters=4000"
    STANDARD_PASS_OUTPUT 
    )

  PACKAGE_COPY_FILES_TO_BINARY_DIR(CopyTestRCGFiles
    SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/epetra/example/RCG
    SOURCE_FILES bcsstk14.hb
    EXEDEPS rcg_hb
    )

ENDIF()
