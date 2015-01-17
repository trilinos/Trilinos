ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LegendreBasisUnitTest
  SOURCES Stokhos_LegendreBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NormalizedLegendreBasisUnitTest
  SOURCES Stokhos_NormalizedLegendreBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  HermiteBasisUnitTest
  SOURCES Stokhos_HermiteBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  NormalizedHermiteBasisUnitTest
  SOURCES Stokhos_NormalizedHermiteBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  JacobiBasisUnitTest
  SOURCES Stokhos_JacobiBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#  DerivExpansionUnitTest
#  SOURCES Stokhos_DerivExpansionUnitTest.cpp
#  COMM serial mpi
#  STANDARD_PASS_OUTPUT
#  NUM_MPI_PROCS 1
#  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  QuadExpansionUnitTest
  SOURCES Stokhos_QuadExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  QuadraturePseudoSpectralExpansionUnitTest
  SOURCES Stokhos_QuadraturePseudoSpectralExpansionUnitTest.cpp
          Stokhos_PseudoSpectralExpansionUnitTest.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorProductPseudoSpectralExpansionUnitTest
  SOURCES Stokhos_TensorProductPseudoSpectralExpansionUnitTest.cpp
          Stokhos_PseudoSpectralExpansionUnitTest.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SmolyakPseudoSpectralExpansionUnitTest
  SOURCES Stokhos_SmolyakPseudoSpectralExpansionUnitTest.cpp
          Stokhos_PseudoSpectralExpansionUnitTest.hpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  AlgebraicExpansionUnitTest
  SOURCES Stokhos_AlgebraicExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ConstantExpansionUnitTest
  SOURCES Stokhos_ConstantExpansionUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  DivisionOperatorUnitTest
  SOURCES Stokhos_DivisionOperatorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

IF(TPL_ENABLE_ForUQTK)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ForUQTKExpansionUnitTest
    SOURCES Stokhos_ForUQTKExpansionUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

ENDIF()

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  StieltjesUnitTest
  SOURCES Stokhos_StieltjesUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LanczosUnitTest
  SOURCES Stokhos_LanczosUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  GramSchmidtUnitTest
  SOURCES Stokhos_GramSchmidtUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Sparse3TensorUnitTest
  SOURCES Stokhos_Sparse3TensorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ExponentialRandomFieldUnitTest
  SOURCES Stokhos_ExponentialRandomFieldUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LogNormalUnitTest
  SOURCES Stokhos_LogNormalUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SDMUtilsUnitTest
  SOURCES Stokhos_SDMUtilsUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  ProductBasisUtilsUnitTest
  SOURCES Stokhos_ProductBasisUtilsUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorProductBasisUnitTest
  SOURCES Stokhos_TensorProductBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TotalOrderBasisUnitTest
  SOURCES Stokhos_TotalOrderBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SmolyakBasisUnitTest
  SOURCES Stokhos_SmolyakBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  TensorProductPseudoSpectralOperatorUnitTest
  SOURCES Stokhos_TensorProductPseudoSpectralOperatorUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  LexicographicTreeBasisUnitTest
  SOURCES Stokhos_LexicographicTreeBasisUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  SparseGridQuadratureUnitTest
  SOURCES Stokhos_SparseGridQuadratureUnitTest.cpp
  COMM serial mpi
  STANDARD_PASS_OUTPUT
  NUM_MPI_PROCS 1
  )

IF(Stokhos_ENABLE_TriKota)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SmolyakPseudoSpectralOperatorUnitTest
    SOURCES Stokhos_SmolyakPseudoSpectralOperatorUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

ENDIF()

IF (Stokhos_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    MatrixFreeOperatorUnitTest
    SOURCES Stokhos_MatrixFreeOperatorUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    InterlacedOpUnitTest
    SOURCES Stokhos_InterlacedOpUnitTest.cpp
            Stokhos_InterlacedTestSupport.cpp
            Stokhos_InterlacedTestSupport.hpp
            ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BasisInteractionGraphUnitTest
    SOURCES Stokhos_BasisInteractionGraphUnitTest.cpp
            Stokhos_InterlacedTestSupport.cpp
            Stokhos_InterlacedTestSupport.hpp
            ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    AdaptivityToolsUnitTest
    SOURCES Stokhos_AdaptivityToolsUnitTest.cpp
            Stokhos_InterlacedTestSupport.cpp
            Stokhos_InterlacedTestSupport.hpp
            ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    BuildColBasisUnitTest
    SOURCES Stokhos_BuildColBasisUnitTest.cpp
            Stokhos_InterlacedTestSupport.cpp
            Stokhos_InterlacedTestSupport.hpp
            ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    InterlacedMapUnitTest
    SOURCES Stokhos_InterlacedMapUnitTest.cpp
            Stokhos_InterlacedTestSupport.cpp
            Stokhos_InterlacedTestSupport.hpp
            ${TEUCHOS_STD_UNIT_TEST_MAIN}
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 2
    )

ENDIF()

IF(Stokhos_ENABLE_Sacado)

  TRIBITS_INCLUDE_DIRECTORIES(REQUIRED_DURING_INSTALLATION_TESTING ${PACKAGE_SOURCE_DIR}/../sacado/test/utils)

  SET(UTILS_SOURCES
    ${PACKAGE_SOURCE_DIR}/../sacado/test/utils/Sacado_Random.hpp
    ${PACKAGE_SOURCE_DIR}/../sacado/test/utils/Sacado_RandomImp.hpp
  )


  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCEUnitTest
    SOURCES Stokhos_SacadoPCEUnitTest.cpp Stokhos_SacadoPCEUnitTest.hpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoETPCEUnitTest
    SOURCES Stokhos_SacadoETPCEUnitTest.cpp Stokhos_SacadoPCEUnitTest.hpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorUnitTest
    SOURCES Stokhos_SacadoVectorUnitTest.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCESerializationTests
    SOURCES Stokhos_SacadoPCESerializationTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoPCECommTests
    SOURCES Stokhos_SacadoPCECommTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorSerializationTests
    SOURCES Stokhos_SacadoVectorSerializationTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoVectorCommTests
    SOURCES Stokhos_SacadoVectorCommTests.cpp ${UTILS_SOURCES}
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    SacadoTraitsUnitTest
    SOURCES Stokhos_SacadoTraitsUnitTest.cpp
    ARGS
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  IF (Stokhos_ENABLE_CXX11)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      SacadoPromoteUnitTest
      SOURCES Stokhos_SacadoPromoteUnitTest.cpp
      ARGS
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
      )
  ENDIF()

  IF (Stokhos_ENABLE_KokkosCore)

    IF (Kokkos_ENABLE_Pthread)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        SacadoMPVectorUnitTest
        SOURCES Stokhos_SacadoMPVectorUnitTest.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        SacadoMPVectorSerializationTests
        SOURCES Stokhos_SacadoMPVectorSerializationTests.cpp ${UTILS_SOURCES}
        ARGS
        COMM serial mpi
        NUM_MPI_PROCS 1
        STANDARD_PASS_OUTPUT
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        SacadoMPVectorCommTests
        SOURCES Stokhos_SacadoMPVectorCommTests.cpp ${UTILS_SOURCES}
        ARGS
        COMM serial mpi
        NUM_MPI_PROCS 1
        STANDARD_PASS_OUTPUT
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewMPVectorUnitTest_Threads
        SOURCES Stokhos_KokkosViewMPVectorUnitTest.hpp
                Stokhos_KokkosViewMPVectorUnitTest_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        SacadoUQPCEUnitTest_Threads
        SOURCES Stokhos_SacadoUQPCEUnitTest.hpp
                Stokhos_SacadoUQPCEUnitTest_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewUQPCEUnitTest_Threads
        SOURCES Stokhos_KokkosViewUQPCEUnitTest.hpp
                Stokhos_KokkosViewUQPCEUnitTest_Threads.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )
    ENDIF()

    IF (Kokkos_ENABLE_OpenMP)

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewMPVectorUnitTest_OpenMP
        SOURCES Stokhos_KokkosViewMPVectorUnitTest.hpp
                Stokhos_KokkosViewMPVectorUnitTest_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewUQPCEUnitTest_OpenMP
        SOURCES Stokhos_KokkosViewUQPCEUnitTest.hpp
                Stokhos_KokkosViewUQPCEUnitTest_OpenMP.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        )

    ENDIF()

    IF (Kokkos_ENABLE_Cuda)
      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewMPVectorUnitTest_Cuda
        SOURCES Stokhos_KokkosViewMPVectorUnitTest.hpp
                Stokhos_KokkosViewMPVectorUnitTest_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )

      TRIBITS_ADD_EXECUTABLE_AND_TEST(
        KokkosViewUQPCEUnitTest_Cuda
        SOURCES Stokhos_KokkosViewUQPCEUnitTest.hpp
                 Stokhos_KokkosViewUQPCEUnitTest_Cuda.cpp
        COMM serial mpi
        STANDARD_PASS_OUTPUT
        NUM_MPI_PROCS 1
        RUN_SERIAL
        )
    ENDIF()

    IF (Stokhos_ENABLE_TpetraKernels)

      IF (HAVE_TPETRA_PTHREAD)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixMPVectorUnitTest_Threads
          SOURCES Stokhos_KokkosCrsMatrixMPVectorUnitTest.hpp
                  Stokhos_KokkosCrsMatrixMPVectorUnitTest_Threads.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          )

        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixUQPCEUnitTest_Threads
          SOURCES Stokhos_KokkosCrsMatrixUQPCEUnitTest.hpp
                  Stokhos_KokkosCrsMatrixUQPCEUnitTest_Threads.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          )
      ENDIF()

      IF (HAVE_TPETRA_OPENMP)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixMPVectorUnitTest_OpenMP
          SOURCES Stokhos_KokkosCrsMatrixMPVectorUnitTest.hpp
                  Stokhos_KokkosCrsMatrixMPVectorUnitTest_OpenMP.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          )

        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixUQPCEUnitTest_OpenMP
          SOURCES Stokhos_KokkosCrsMatrixUQPCEUnitTest.hpp
                  Stokhos_KokkosCrsMatrixUQPCEUnitTest_OpenMP.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          )
      ENDIF()

      IF (HAVE_TPETRA_CUDA)
        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixMPVectorUnitTest_Cuda
          SOURCES Stokhos_KokkosCrsMatrixMPVectorUnitTest.hpp
                  Stokhos_KokkosCrsMatrixMPVectorUnitTest_Cuda.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          RUN_SERIAL
          )

        TRIBITS_ADD_EXECUTABLE_AND_TEST(
          KokkosCrsMatrixUQPCEUnitTest_Cuda
          SOURCES Stokhos_KokkosCrsMatrixUQPCEUnitTest.hpp
                  Stokhos_KokkosCrsMatrixUQPCEUnitTest_Cuda.cpp
          COMM serial mpi
          STANDARD_PASS_OUTPUT
          NUM_MPI_PROCS 1
          RUN_SERIAL
          )
      ENDIF()

      IF (Stokhos_ENABLE_Tpetra AND Tpetra_Have_Kokkos_Refactor)

        IF (HAVE_TPETRA_PTHREAD)
          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixMPVectorUnitTest_Threads
            SOURCES Stokhos_TpetraCrsMatrixMPVectorUnitTest.hpp
                    Stokhos_TpetraCrsMatrixMPVectorUnitTest_Threads.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            )

          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixUQPCEUnitTest_Threads
            SOURCES Stokhos_TpetraCrsMatrixUQPCEUnitTest.hpp
                    Stokhos_TpetraCrsMatrixUQPCEUnitTest_Threads.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            )
        ENDIF()

        IF (HAVE_TPETRA_OPENMP)
          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixMPVectorUnitTest_OpenMP
            SOURCES Stokhos_TpetraCrsMatrixMPVectorUnitTest.hpp
                    Stokhos_TpetraCrsMatrixMPVectorUnitTest_OpenMP.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            )

          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixUQPCEUnitTest_OpenMP
            SOURCES Stokhos_TpetraCrsMatrixUQPCEUnitTest.hpp
                    Stokhos_TpetraCrsMatrixUQPCEUnitTest_OpenMP.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            )
        ENDIF()

        IF (HAVE_TPETRA_CUDA)
          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixMPVectorUnitTest_Cuda
            SOURCES Stokhos_TpetraCrsMatrixMPVectorUnitTest.hpp
                    Stokhos_TpetraCrsMatrixMPVectorUnitTest_Cuda.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            RUN_SERIAL
            )

          TRIBITS_ADD_EXECUTABLE_AND_TEST(
            TpetraCrsMatrixUQPCEUnitTest_Cuda
            SOURCES Stokhos_TpetraCrsMatrixUQPCEUnitTest.hpp
                    Stokhos_TpetraCrsMatrixUQPCEUnitTest_Cuda.cpp
            COMM serial mpi
            STANDARD_PASS_OUTPUT
            NUM_MPI_PROCS 4
            RUN_SERIAL
            )
        ENDIF()

        TRIBITS_COPY_FILES_TO_BINARY_DIR(copyMueLuFiles
          DEST_FILES muelu_cheby.xml
          )

      ENDIF()

    ENDIF()

  ENDIF()

ENDIF()

# Need to fix these tests so they don't require Epetra
IF (Stokhos_ENABLE_KokkosCore AND Stokhos_ENABLE_EpetraExt)

#  TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    KokkosArrayKernelsUnitTest_Serial
#    SOURCES Stokhos_KokkosArrayKernelsUnitTest.hpp
#            Stokhos_KokkosArrayKernelsUnitTestDecl.hpp
#            Stokhos_KokkosArrayKernelsUnitTest_Host.hpp
#            Stokhos_KokkosArrayKernelsUnitTest_Serial.cpp
#    COMM serial mpi
#    STANDARD_PASS_OUTPUT
#    NUM_MPI_PROCS 1
##    )

  IF (Kokkos_ENABLE_Pthread)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      KokkosArrayKernelsUnitTest_Threads
      SOURCES Stokhos_KokkosArrayKernelsUnitTest.hpp
              Stokhos_KokkosArrayKernelsUnitTestDecl.hpp
              Stokhos_KokkosArrayKernelsUnitTest_Host.hpp
              Stokhos_KokkosArrayKernelsUnitTest_Threads.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )
  ENDIF()

  IF (Kokkos_ENABLE_OpenMP)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      KokkosArrayKernelsUnitTest_OpenMP
      SOURCES Stokhos_KokkosArrayKernelsUnitTest.hpp
              Stokhos_KokkosArrayKernelsUnitTestDecl.hpp
              Stokhos_KokkosArrayKernelsUnitTest_Host.hpp
              Stokhos_KokkosArrayKernelsUnitTest_OpenMP.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      )
  ENDIF()

  IF (Kokkos_ENABLE_Cuda)
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      KokkosArrayKernelsUnitTest_Cuda
      SOURCES Stokhos_KokkosArrayKernelsUnitTest.hpp
              Stokhos_KokkosArrayKernelsUnitTestDecl.hpp
              Stokhos_KokkosArrayKernelsUnitTest_Cuda.cpp
      COMM serial mpi
      STANDARD_PASS_OUTPUT
      NUM_MPI_PROCS 1
      RUN_SERIAL
      )
  ENDIF()

ENDIF()

IF (Stokhos_ENABLE_Tpetra)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    TpetraMatVec
    SOURCES tpetra_mat_vec.cpp
    COMM serial mpi
    STANDARD_PASS_OUTPUT
    )
ENDIF()
