TRIBITS_INCLUDE_DIRECTORIES(${Belos_SOURCE_DIR}/kokkos/src)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_BlockGmresEx
  SOURCES BlockGmresKokkosExFile.cpp
  ARGS "--verbose"
       "--verbose --num-rhs=4"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_KKGMRESEx
  SOURCES ex_real_A.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_BlockCGEx
  SOURCES BlockCGKokkosExFile.cpp
  ARGS "--verbose"
       "--verbose --num-rhs=4"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_GCRODREx
  SOURCES GCRODRKokkosExFile.cpp
  ARGS "--verbose"
  COMM serial mpi
  NUM_MPI_PROCS 1
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleMtxFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
  SOURCE_FILES bcsstk12.mtx
  EXEDEPS Kokkos_BlockCGEx Kokkos_BlockGmresEx Kokkos_GCRODREx
  )

IF(Teuchos_ENABLE_COMPLEX)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    Kokkos_GCRODRComplexEx
    SOURCES GCRODRComplexKokkosExFile.cpp
    ARGS "--verbose"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleComplexMtxFiles
    SOURCE_DIR ${Belos_SOURCE_DIR}/testmatrices
    SOURCE_FILES mhd1280b.mtx
    EXEDEPS Kokkos_GCRODRComplexEx
    )

ENDIF(Teuchos_ENABLE_COMPLEX)
