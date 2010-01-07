
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

# PACKAGE_ADD_EXECUTABLE_AND_TEST(
#   BasicPerfTest
#   SOURCES
#     BasicPerfTest
#     ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
#   # ARGS 
#   COMM serial # mpi
#   STANDARD_PASS_OUTPUT
#   )


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  HybridPlatformTest
  SOURCES
    HybridPlatformTest
  # ARGS 
  COMM mpi # serial
  STANDARD_PASS_OUTPUT
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files1
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/tpetra/test/BasicPerfTest
  SOURCE_FILES machine.xml
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(Tpetra_HybridPlatformTest_Files2
  SOURCE_DIR ${CMAKE_SOURCE_DIR}/packages/belos/tpetra/test/BlockCG
  SOURCE_FILES bcsstk17.rsa
)
