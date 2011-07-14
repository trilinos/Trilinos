INCLUDE(PackageAddExecutableAndTest)


PACKAGE_ADD_EXECUTABLE_AND_TEST(
  NodeAPI_DefaultNodeTest
  SOURCES 
    DefaultNodeTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
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
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkosnodeapi
  )
SET(TIMING_INSTALLS Kokkos_SerialNodeTestAndTiming)

IF (Kokkos_ENABLE_TBB)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TBBNodeTestAndTiming
    SOURCES 
      NodeTest
      TBBNodeTest
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
  APPEND_SET(TIMING_INSTALLS Kokkos_TBBNodeTestAndTiming)
ENDIF()

IF (Kokkos_ENABLE_ThreadPool)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    TPINodeTestAndTiming
    SOURCES 
      NodeTest
      TPINodeTest
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
  APPEND_SET(TIMING_INSTALLS Kokkos_TPINodeTestAndTiming)
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
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
  APPEND_SET(TIMING_INSTALLS Kokkos_ThrustGPUNodeTestAndTiming)

  CUDA_ADD_EXECUTABLE(Kokkos_CUDABuiltAxpy CUDABuiltAxpy.cu)
  TARGET_LINK_LIBRARIES(Kokkos_CUDABuiltAxpy kokkosnodeapi)
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
