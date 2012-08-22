INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsAddExecutable)

# As of <= 14 Mar 2012, MultiVectorFiller doesn't work yet with
# multiple MPI processes.  Thus, I'll leave it out of the MPI tests
# for now.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MultiVectorFiller_SerialNodeTest
  SOURCES MultiVectorFiller_SerialNodeTest
  COMM serial
  STANDARD_PASS_OUTPUT
  DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  )

IF (KokkosClassic_ENABLE_TBB AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_TbbTest
  #   SOURCES MultiVectorFiller_TbbTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

IF (KokkosClassic_ENABLE_ThreadPool AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThreadPoolTest
  #   SOURCES MultiVectorFiller_ThreadPoolTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

IF (KokkosClassic_ENABLE_Thrust AND KokkosClassic_ENABLE_CUDA_DOUBLE AND Tpetra_ENABLE_MPI)
  # TRIBITS_ADD_EXECUTABLE_AND_TEST(
  #   MultiVectorFiller_ThrustTest
  #   SOURCES MultiVectorFiller_ThrustTest
  #   COMM mpi
  #   STANDARD_PASS_OUTPUT
  #   DEPLIBS kokkos kokkoslinalg kokkosnodeapi
  #   )
ENDIF()

