INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# IF (Anasazi_ENABLE_Tpetra)
# ENDIF()

# IF (Belos_ENABLE_Tpetra)
# ENDIF()


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
    HybridPlatformTest
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi tpetrainout tpetratestlib
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/tpetra/example/HybridPlatform
  SOURCE_FILES gpuonly.xml hybrid.xml mpionly.xml threaded.xml
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)
