TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/test/utils)

SET(UTILS_SOURCES
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_Random.hpp
  ${PACKAGE_SOURCE_DIR}/test/utils/Sacado_RandomImp.hpp
)

# All of these tests use the Teuchos test harness, and therefore TeuchosCore
# must be enabled
IF (Sacado_ENABLE_TeuchosCore)

  IF (Sacado_ENABLE_TeuchosComm)
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

  IF (Sacado_ENABLE_Kokkos)

    IF (Kokkos_ENABLE_SERIAL)
      IF (Sacado_ENABLE_TeuchosKokkosComm)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Fad_Kokkos_CommTests_Serial
          SOURCES Fad_Kokkos_CommTests_Serial.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
          ARGS
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          )
      ENDIF()

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Serial
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Serial.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )
      
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
	FadKokkosTests_ViewCornerCases
	SOURCES Fad_KokkosTests_ViewCornerCases.cpp
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

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadFadKokkosTests_Serial
        SOURCES Fad_Fad_KokkosTests.hpp
                Fad_Fad_KokkosTests_Serial.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Serial
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Serial.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_NoViewSpec_Serial
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_NoViewSpec_Serial.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )
    ENDIF ()

    IF (Kokkos_ENABLE_THREADS)
      IF (Sacado_ENABLE_TeuchosKokkosComm)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Fad_Kokkos_CommTests_Threads
          SOURCES Fad_Kokkos_CommTests_Threads.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
          ARGS
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          )
      ENDIF()

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Threads
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_NoViewSpec_Threads
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_NoViewSpec_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadFadKokkosTests_Threads
        SOURCES Fad_Fad_KokkosTests.hpp
                Fad_Fad_KokkosTests_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Threads
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_NoViewSpec_Threads
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_NoViewSpec_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )
    ENDIF()

    IF (Kokkos_ENABLE_OPENMP)
      IF (Sacado_ENABLE_TeuchosKokkosComm)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Fad_Kokkos_CommTests_OpenMP
          SOURCES Fad_Kokkos_CommTests_OpenMP.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
          ARGS
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          )
      ENDIF()

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_OpenMP
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_NoViewSpec_OpenMP
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_NoViewSpec_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadFadKokkosTests_OpenMP
        SOURCES Fad_Fad_KokkosTests.hpp
                Fad_Fad_KokkosTests_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_OpenMP
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_NoViewSpec_OpenMP
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_NoViewSpec_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )
    ENDIF()

    IF (Kokkos_ENABLE_CUDA)
      IF (Sacado_ENABLE_TeuchosKokkosComm)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Fad_Kokkos_CommTests_Cuda
          SOURCES Fad_Kokkos_CommTests_Cuda.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
          ARGS
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          )
      ENDIF()

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Cuda
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Cuda_Hierarchical
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Cuda_Hierarchical.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Cuda_Hierarchical_SFad
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Cuda_Hierarchical_SFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Cuda_Hierarchical_DFad
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_Cuda_Hierarchical_DFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_NoViewSpec_Cuda
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_NoViewSpec_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadFadKokkosTests_Cuda
        SOURCES Fad_Fad_KokkosTests.hpp
                Fad_Fad_KokkosTests_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Cuda
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_NoViewSpec_Cuda
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_NoViewSpec_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Cuda_Hierarchical
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Cuda_Hierarchical.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Cuda_Hierarchical_SFad
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Cuda_Hierarchical_SFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_Cuda_Hierarchical_DFad
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_Cuda_Hierarchical_DFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )
    ENDIF()

    IF (Kokkos_ENABLE_HIP)
      IF (Sacado_ENABLE_TeuchosKokkosComm)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          Fad_Kokkos_CommTests_HIP
          SOURCES Fad_Kokkos_CommTests_HIP.cpp Fad_CommTests.hpp ${UTILS_SOURCES}
          ARGS
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          )
      ENDIF()

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_HIP
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_HIP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_HIP_Hierarchical
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_HIP_Hierarchical.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_HIP_Hierarchical_SFad
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_HIP_Hierarchical_SFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_Cuda_Hierarchical_DFad
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_HIP_Hierarchical_DFad.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosTests_NoViewSpec_HIP
        SOURCES Fad_KokkosTests.hpp
                Fad_KokkosTests_NoViewSpec_HIP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        FadFadKokkosTests_Cuda
#        SOURCES Fad_Fad_KokkosTests.hpp
#                Fad_Fad_KokkosTests_Cuda.cpp
#        COMM serial mpi
#        STANDARD_PASS_OUTPUT
#        NUM_MPI_PROCS 1
#        RUN_SERIAL
#        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_HIP
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_HIP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        FadKokkosAtomicTests_NoViewSpec_HIP
        SOURCES Fad_KokkosAtomicTests.hpp
                Fad_KokkosAtomicTests_NoViewSpec_HIP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        FadKokkosAtomicTests_Cuda_Hierarchical
#        SOURCES Fad_KokkosAtomicTests.hpp
#                Fad_KokkosAtomicTests_Cuda_Hierarchical.cpp
#        COMM serial mpi
#        STANDARD_PASS_OUTPUT
#        NUM_MPI_PROCS 1
#        RUN_SERIAL
#        )

#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        FadKokkosAtomicTests_Cuda_Hierarchical_SFad
#        SOURCES Fad_KokkosAtomicTests.hpp
#                Fad_KokkosAtomicTests_Cuda_Hierarchical_SFad.cpp
#        COMM serial mpi
#        STANDARD_PASS_OUTPUT
#        NUM_MPI_PROCS 1
#        RUN_SERIAL
#        )

#      TRIBITS_ADD_EXECUTABLE_AND_TEST(
#        FadKokkosAtomicTests_Cuda_Hierarchical_DFad
#        SOURCES Fad_KokkosAtomicTests.hpp
#                Fad_KokkosAtomicTests_Cuda_Hierarchical_DFad.cpp
#        COMM serial mpi
#        STANDARD_PASS_OUTPUT
#        NUM_MPI_PROCS 1
#        RUN_SERIAL
#        )
    ENDIF()

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      ViewFactoryTests
      SOURCES ViewFactoryTests.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      RUN_SERIAL
      )

  ENDIF()

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    PromoteTests
    SOURCES PromoteTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MPLTests
    SOURCES MPLTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ConditionalReturnTypeTest
    SOURCES ConditionalReturnTypeTest.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ConversionTests
    SOURCES ConversionTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MoveTests
    SOURCES MoveTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TernaryExprTests
    SOURCES TernaryExprTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SafeSqrtTests
    SOURCES SafeSqrtTests.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    )
ENDIF()
