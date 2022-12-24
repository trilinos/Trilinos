INCLUDE_DIRECTORIES(${Belos_SOURCE_DIR}/kokkos/src)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_BlockGmresEx
  SOURCES BlockGmresKokkosExFile.cpp
  ARGS "--verbose"
       "--verbose --num-rhs=4"
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  Kokkos_BlockCGEx
  SOURCES BlockCGKokkosExFile.cpp
  ARGS "--verbose"
       "--verbose --num-rhs=4"
  COMM serial mpi
  )

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyExampleMtxFiles
  SOURCE_DIR ${Belos_SOURCE_DIR}/kokkos/example
  SOURCE_FILES bcsstk12.mtx
  EXEDEPS Kokkos_BlockCGEx Kokkos_BlockGmresEx
  )
