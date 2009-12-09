
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/NodeAPI)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../kokkos/LinAlg)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../kokkos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../tpetra/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../tpetra/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../belos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../belos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../belos/tpetra/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../belos/tpetra/src)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

SET(SOURCES "")

APPEND_SET(SOURCES
  belos_solve.cpp
)

PACKAGE_ADD_EXECUTABLE(
  Tifpack_belos
  SOURCES ${SOURCES}
  DEPLIBS belos tpetra kokkos teuchos
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(TifpackBelosCopyFiles
  DEST_FILES
  ILUT_calore1_mm.xml PointJacobi_calore1_mm.xml GaussSeidel_calore1_mm.xml
  SymGaussSeidel_calore1_mm.xml calore1.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS Tifpack_belos
)

PACKAGE_ADD_TEST(
  Tifpack_belos TifpackBelosCopyFiles
  NAME PointJacobi_belos
  ARGS "--xml_params_file=PointJacobi_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  Tifpack_belos TifpackBelosCopyFiles
  NAME GaussSeidel_belos
  ARGS "--xml_params_file=GaussSeidel_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  Tifpack_belos TifpackBelosCopyFiles
  NAME SymGaussSeidel_belos
  ARGS "--xml_params_file=SymGaussSeidel_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  Tifpack_belos TifpackBelosCopyFiles
  NAME ILUT_belos
  ARGS "--xml_params_file=ILUT_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

