
INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)


TRIBITS_ADD_EXECUTABLE(
  ifp_aztecoo
  SOURCES aztecoo_solve.cpp build_problem.cpp build_solver.cpp read_matrix.cpp
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(IfpackAztecOOCopyFiles
  DEST_FILES
  test_ILUT_tfqmr_small_sym_mm.xml small_sym.mtx test_gmres_small_sym_mm.xml
  test_az_ILUT_tfqmr_small_sym_mm.xml
  test_2_ILUT_nos1_hb.xml nos1.rsa
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS ifp_aztecoo
)

TRIBITS_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_gmres_az
  ARGS "--xml_file=test_gmres_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_ILUT_tfqmr
  ARGS "--xml_file=test_ILUT_tfqmr_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_az_ILUT_tfqmr
  ARGS "--xml_file=test_az_ILUT_tfqmr_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME small_az_ILUT_tfqmr_p3
  ARGS "--xml_file=test_az_ILUT_tfqmr_small_sym_mm.xml"
  COMM mpi
  NUM_MPI_PROCS 3
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  ifp_aztecoo IfpackAztecOOCopyFiles
  NAME ILUT_nos1_hb
  ARGS "--xml_file=test_2_ILUT_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

