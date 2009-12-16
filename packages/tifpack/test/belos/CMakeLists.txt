
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


PACKAGE_ADD_EXECUTABLE(
  tif_belos
  SOURCES belos_solve.cpp
  DEPLIBS belos tpetra kokkos teuchos
)

IF(Tifpack_ENABLE_QD)
  PACKAGE_ADD_EXECUTABLE(
    tif_extprec_belos
    SOURCES belos_extprec_solve.cpp
    DEPLIBS belos tpetra kokkos teuchos
  )
ENDIF()

PACKAGE_COPY_FILES_TO_BINARY_DIR(TifpackBelosCopyFiles
  DEST_FILES
  test_ILUT_calore1_mm.xml test_PtJacobi_calore1_mm.xml test_GS_calore1_mm.xml
  test_SGS_calore1_mm.xml calore1.mtx
  test_PtJacobi_small_sym_mm.xml small_sym.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tif_belos
)

PACKAGE_ADD_TEST(
  tif_belos TifpackBelosCopyFiles
  NAME PtJacobi_belos
  ARGS "--xml_params_file=test_PtJacobi_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

# IF(Tifpack_ENABLE_QD)
#   PACKAGE_ADD_TEST(
#     tif_extprec_belos TifpackBelosCopyFiles
#     NAME PtJacobi_extprec_belos
#     ARGS "--xml_params_file=test_PtJacobi_small_sym_mm.xml"
#     COMM serial mpi
#     NUM_MPI_PROCS 1
#     STANDARD_PASS_OUTPUT
#   )
# ENDIF()

PACKAGE_ADD_TEST(
  tif_belos TifpackBelosCopyFiles
  NAME PtJacobi_belos
  ARGS "--xml_params_file=test_PtJacobi_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  tif_belos TifpackBelosCopyFiles
  NAME GS_belos
  ARGS "--xml_params_file=test_GS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  tif_belos TifpackBelosCopyFiles
  NAME SGS_belos
  ARGS "--xml_params_file=test_SGS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  tif_belos TifpackBelosCopyFiles
  NAME ILUT_belos
  ARGS "--xml_params_file=test_ILUT_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

