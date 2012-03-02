INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddExecutable)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVectorFiller_SerialNodeTest
  SOURCES MultiVectorFiller_SerialNodeTest
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

IF (Kokkos_ENABLE_TBB AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_TbbTest
  #   SOURCES MultiVectorFiller_TbbTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

IF (Kokkos_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThreadPoolTest
  #   SOURCES MultiVectorFiller_ThreadPoolTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

IF (Kokkos_ENABLE_Thrust AND Kokkos_ENABLE_CUDA_DOUBLE AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThrustTest
  #   SOURCES MultiVectorFiller_ThrustTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

