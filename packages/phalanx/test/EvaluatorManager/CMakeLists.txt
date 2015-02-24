

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR})

TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${CMAKE_CURRENT_SOURCE_DIR}/../Utilities)

IF(PHX_KOKKOS_DEVICE_TYPE_CUDA2)
  IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
  ENDIF()

  CUDA_INCLUDE_DIRECTORIES( ${KokkosCore_INCLUDE_DIRS}
                            ${CMAKE_CURRENT_BINARY_DIR}
                            ${CMAKE_CURRENT_SOURCE_DIR}/../../src 
                            ${CMAKE_CURRENT_BINARY_DIR}/../../src
                            ${CMAKE_CURRENT_SOURCE_DIR}/../../../sacado/src
                            ${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/core/src
                            ${Phalanx_INCLUDE_DIRS}
                            )
  CUDA_COMPILE( OBJECTS_CUDA EvaluatorManager_CU.cu )

  LIST( APPEND SOURCES ${OBJECTS_CUDA} )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    evaluator_manager
    SOURCES ${SOURCES} ${TEUCHOS_STD_UNIT_TEST_MAIN}
    TESTONLYLIBS phalanx_test_utilities
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

ELSE()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    evaluator_manager
    SOURCES EvaluatorManager.cpp ${TEUCHOS_STD_UNIT_TEST_MAIN}
    TESTONLYLIBS phalanx_test_utilities
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

ENDIF ()

