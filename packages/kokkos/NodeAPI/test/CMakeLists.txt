INCLUDE(PackageAddExecutableAndTest)

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
  SerialNodeTestAndTiming
  SOURCES 
    NodeTest
    SerialNodeTest
    ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkosnodeapi
  )

IF (Kokkos_ENABLE_TBB)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TBBNodeTestAndTiming
    SOURCES 
      NodeTest
      TBBNodeTest
      ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
ENDIF()

IF (Kokkos_ENABLE_ThreadPool)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TPINodeTestAndTiming
    SOURCES 
      NodeTest
      TPINodeTest
      ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
ENDIF()

IF (Kokkos_ENABLE_Thrust)
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(source_files TestOps.cu)
  CUDA_COMPILE(CUDA_FILES TestOps.cu)

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    ThrustGPUNodeTestAndTiming
    SOURCES 
      NodeTest
      ThrustGPUNodeTest
      ${CUDA_FILES}
      ${source_files}
      ${PACKAGE_SOURCE_DIR}/../teuchos/test/UnitTest/Teuchos_StandardUnitTestMain.cpp
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
ENDIF()
