INCLUDE(PackageAddExecutableAndTest)

ASSERT_DEFINED(Kokkos_ENABLE_Thrust)
IF (Kokkos_ENABLE_Thrust)
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(source_files TestOps.cu)
  CUDA_COMPILE(CUDA_FILES TestOps.cu)
ENDIF()

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NodeAPI_DefaultNodeTest
  SOURCES 
    DefaultNodeTest
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkosnodeapi
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NodeAPI_UnitTests
  SOURCES 
    NodeTest
    ${CUDA_FILES}
    ${source_files}
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkosnodeapi
  )
