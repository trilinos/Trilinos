INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

IF (Trilinos_ENABLE_CXX11)
  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    RTIExample
    SOURCES RTIExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    RTIOperatorExample
    SOURCES RTIOperatorExample
    ARGS 
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  PACKAGE_ADD_EXECUTABLE_AND_TEST(
    RTIInlineCG
    SOURCES RTIInlineCG
    ARGS 
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  IF (Tpetra_ENABLE_QD)
    PACKAGE_ADD_EXECUTABLE_AND_TEST(
      RTIInlineCG_qd
      SOURCES RTIInlineCG_qd
      ARGS 
      COMM serial mpi
      STANDARD_PASS_OUTPUT
    )
  ENDIF()

# ASSERT_DEFINED(Kokkos_ENABLE_Thrust)
# IF (Kokkos_ENABLE_Thrust)
#   # this pulls in Thrust, CUSP, and Kokkos/NodeAPI
#   include_directories(${Kokkos_INCLUDE_DIRS})
#   CUDA_ADD_EXECUTABLE(Tpetra_RTIExampleThrust RTIExampleThrust.cu)
#   PACKAGE_ADD_TEST(
#     Tpetra_RTIExampleThrust 
#     ARGS
#     COMM serial mpi
#     NUM_MPI_PROCS 1
#     STANDARD_PASS_OUTPUT
#   )
# ENDIF()

ENDIF()
