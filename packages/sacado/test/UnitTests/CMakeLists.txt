TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/utils)

SET(UTILS_SOURCES
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_Random.hpp
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_RandomImp.hpp
)

IF (Sacado_ENABLE_Teuchos)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FadSerializationTests
    SOURCES Fad_SerializationTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TaySerializationTests
    SOURCES Tay_SerializationTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    XHOSTTYPE Windows Windows_NT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    FadCommTests
    SOURCES Fad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ELRFadCommTests
    SOURCES ELRFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    CacheFadCommTests
    SOURCES CacheFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ELRCacheFadCommTests
    SOURCES ELRCacheFad_CommTests.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TayCommTests
    SOURCES Tay_CommTests.cpp Tay_CommTests.hpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    XHOSTTYPE Windows Windows_NT
  )
ENDIF()

IF (Sacado_ENABLE_KokkosCore)

  IF (Kokkos_ENABLE_Serial)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_Serial
      SOURCES Fad_KokkosTests.hpp
        Fad_KokkosTests_Serial.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_NoViewSpec_Serial
      SOURCES Fad_KokkosTests.hpp
        Fad_KokkosTests_NoViewSpec_Serial.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )
  ENDIF ()

  IF (Kokkos_ENABLE_Pthread)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_Threads
      SOURCES Fad_KokkosTests.hpp
              Fad_KokkosTests_Threads.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_NoViewSpec_Threads
      SOURCES Fad_KokkosTests.hpp
              Fad_KokkosTests_NoViewSpec_Threads.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )
  ENDIF()

  IF (Kokkos_ENABLE_OpenMP)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_OpenMP
      SOURCES Fad_KokkosTests.hpp
              Fad_KokkosTests_OpenMP.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_NoViewSpec_OpenMP
      SOURCES Fad_KokkosTests.hpp
              Fad_KokkosTests_NoViewSpec_OpenMP.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )
  ENDIF()

  IF (Kokkos_ENABLE_Cuda)
    IF( TPL_ENABLE_MPI )
      CUDA_INCLUDE_DIRECTORIES( ${MPI_BASE_DIR}/include )
    ENDIF()
    CUDA_COMPILE(CUDA_SOURCES Fad_KokkosTests_Cuda.cu)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_Cuda
      SOURCES Fad_KokkosTests.hpp
              ${CUDA_SOURCES}
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      RUN_SERIAL
      )

    CUDA_COMPILE(CUDA_SOURCES Fad_KokkosTests_NoViewSpec_Cuda.cu)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      FadKokkosTests_NoViewSpec_Cuda
      SOURCES Fad_KokkosTests.hpp
              ${CUDA_SOURCES}
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      RUN_SERIAL
      )
  ENDIF()

ENDIF()

IF (Sacado_ENABLE_CXX11)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    PromoteTests
    SOURCES PromoteTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConversionTests
  SOURCES ConversionTests.cpp
  ARGS
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )
