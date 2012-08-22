INCLUDE(TribitsAddExecutableAndTest)

INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NodeAPI_DefaultNodeTest
  SOURCES 
    DefaultNodeTest
    ${TEUCHOS_STD_UNIT_TEST_MAIN}
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkosnodeapi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
SET(TIMING_INSTALLS KokkosClassic_SerialNodeTestAndTiming)

IF (KokkosClassic_ENABLE_TBB)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
  APPEND_SET(TIMING_INSTALLS KokkosClassic_TBBNodeTestAndTiming)
ENDIF()

IF (KokkosClassic_ENABLE_OpenMP)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    OpenMPNodeTestAndTiming
    SOURCES 
      NodeTest
      OpenMPNodeTest
      ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    DEPLIBS kokkosnodeapi
    )
  APPEND_SET(TIMING_INSTALLS KokkosClassic_OpenMPNodeTestAndTiming)
ENDIF()

IF (KokkosClassic_ENABLE_ThreadPool)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
  APPEND_SET(TIMING_INSTALLS KokkosClassic_TPINodeTestAndTiming)
ENDIF()

IF (KokkosClassic_ENABLE_Thrust)
  set(CUDA_ATTACH_VS_BUILD_RULE_TO_CUDA_FILE OFF)
  set(source_files TestOps.cu)
  CUDA_COMPILE(CUDA_FILES TestOps.cu)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
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
  APPEND_SET(TIMING_INSTALLS KokkosClassic_ThrustGPUNodeTestAndTiming)

  CUDA_COMPILE(CUDA_FILES CUDABuiltAxpy.cu)
  TRIBITS_ADD_EXECUTABLE(
    CUDABuiltAxpy 
    SOURCES CUDABuiltAxpy.cpp
            ${CUDA_FILES}
    COMM serial mpi
    DEPLIBS kokkosnodeapi
  )
ENDIF()

INSTALL(TARGETS ${TIMING_INSTALLS}
        RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}")
