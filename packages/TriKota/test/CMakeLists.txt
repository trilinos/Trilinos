
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/install/lib/)


# This is needed to compile codes with dakota headers.
ADD_DEFINITIONS(-DHAVE_CONFIG_H -DDISABLE_DAKOTA_CONFIG_H)

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleEpetraExtME
  SOURCES 
  Main_SimpleEpetraExtME.cpp
  Simple_EpetraExtME.cpp
  Simple_EpetraExtME.hpp
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TriKotaSimpleCopyDakotaIn
  DEST_FILES   dakota.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS  SimpleEpetraExtME
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  SimpleThyraME
  SOURCES 
  Main_SimpleThyraME.cpp
  Simple_EpetraExtME.cpp
  Simple_EpetraExtME.hpp
  COMM serial mpi
  NUM_MPI_PROCS 2
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TriKotaThyraMECopyDakotaIn
  DEST_FILES   dakota.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS  SimpleThyraME
  )

PACKAGE_ADD_EXECUTABLE_AND_TEST(
  ParallelDiagonalThyraME 
  SOURCES 
  Main_ParallelDiagonalThyraME.cpp
  Diagonal_ThyraROME_def.hpp
  Diagonal_ThyraROME.hpp
  COMM serial mpi
  NUM_MPI_PROCS 4
  ARGS -v
  PASS_REGULAR_EXPRESSION "TEST PASSED"
  )

PACKAGE_COPY_FILES_TO_BINARY_DIR(TriKotaParallelDiagonalThyraMECopyDakotaIn
  DEST_FILES   dakota_conmin.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS ParallelDiagonalThyraME
  )
