
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
LINK_DIRECTORIES(${PACKAGE_BINARY_DIR}/Dakota/install/lib/)

# Needed to find Dakota components built with ExternalProject_Add
# and their dependencies; only defined in recent versions of Dakota
IF(DEFINED Dakota_LINK_DIRS)
  LINK_DIRECTORIES(${Dakota_LINK_DIRS})
ENDIF()

# This is needed to compile codes with dakota headers.
ADD_DEFINITIONS(${TRIKOTA_DAKOTA_DEFINITIONS})

TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

TRIBITS_COPY_FILES_TO_BINARY_DIR(TriKotaSimpleCopyDakotaIn
  DEST_FILES   dakota.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS  SimpleEpetraExtME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

TRIBITS_COPY_FILES_TO_BINARY_DIR(TriKotaThyraMECopyDakotaIn
  DEST_FILES   dakota.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS  SimpleThyraME
  )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
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

TRIBITS_COPY_FILES_TO_BINARY_DIR(TriKotaParallelDiagonalThyraMECopyDakotaIn
  DEST_FILES   dakota_conmin.in
  SOURCE_DIR   ${PACKAGE_SOURCE_DIR}/test
  SOURCE_PREFIX "_"
  EXEDEPS ParallelDiagonalThyraME
  )
