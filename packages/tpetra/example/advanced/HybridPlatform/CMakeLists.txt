TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
  HybridPlatformTest
  # ARGS
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_EXECUTABLE(
  HybridPlatformTutorial
  SOURCES HybridPlatformTutorial
  COMM mpi # serial
  )

# mfh 19 Oct 2014: This test only runs correctly if all the given Node
# types are enabled.  However, the test should still build, which is
# why we separate the executable declaration above from the test
# declaration below.
IF (KokkosClassic_ENABLE_Serial AND KokkosClassic_ENABLE_OpenMP AND KokkosClassic_ENABLE_TBB AND KokkosClassic_ENABLE_ThreadPool AND Kokkos_ENABLE_CUDA AND KokkosClassic_ENABLE_CUDA_DOUBLE)
  TRIBITS_ADD_TEST(
    HybridPlatformTutorial
    NAME HybridPlatformTutorial
    # ARGS
    COMM mpi # serial
    STANDARD_PASS_OUTPUT
    COMM mpi
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${Tpetra_MACHINE_XML_FILE_DIR}
  SOURCE_FILES default.xml gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${Belos_SOURCE_DIR}/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)
