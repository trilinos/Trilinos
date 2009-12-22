
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../teuchos/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../aztecoo/src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/../../../epetra/src)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


PACKAGE_ADD_EXECUTABLE(
  ifp_aztecoo
  SOURCES aztecoo_solve.cpp build_problem.cpp build_solver.cpp read_matrix.cpp
  DEPLIBS aztecoo epetra teuchos
)

PACKAGE_COPY_FILES_TO_BINARY_DIR(IfpackAztecOOCopyFiles
  DEST_FILES
  test_ILUT_tfqmr_small_sym_mm.xml small_sym.mtx test_gmres_small_sym_mm.xml
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS ifp_aztecoo
)

PACKAGE_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_gmres_az
  ARGS "--xml_params_file=test_gmres_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

PACKAGE_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_ILUT_tfqmr
  ARGS "--xml_params_file=test_ILUT_tfqmr_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

