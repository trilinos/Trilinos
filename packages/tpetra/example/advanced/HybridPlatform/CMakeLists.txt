
# These tests require double, so disable if cuda is enabled
# but cuda-double is not enabled
SET(do_hybrid_tests ON)
IF( Kokkos_ENABLE_Cuda AND NOT KokkosClassic_ENABLE_CUDA_DOUBLE )
  SET(do_hybrid_tests OFF)
ENDIF()
IF( do_hybrid_tests )
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    HybridPlatformTest
    SOURCES
      HybridPlatformTest
    # ARGS
    COMM mpi # serial
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    HybridPlatformTutorial
    SOURCES
      HybridPlatformTutorial
    # ARGS
    COMM mpi # serial
    STANDARD_PASS_OUTPUT
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${Tpetra_MACHINE_XML_FILE_DIR}
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)
