INCLUDE(TribitsAddExecutableAndTest)
INCLUDE(TribitsCopyFilesToBinaryDir)

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_SUBDIRECTORY(hybrid_platform)

TRIBITS_ADD_EXECUTABLE(
  tif_belos
  SOURCES belos_solve.cpp
)

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_EXECUTABLE(
    tif_extprec_belos
    SOURCES belos_extprec_solve.cpp
    DEPLIBS ifpack2_tpetra_testinst
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2BelosCopyFiles
  DEST_FILES
  test_ILUT_calore1_mm.xml test_Jacobi_calore1_mm.xml test_GS_calore1_mm.xml
  test_RILUK_calore1_mm.xml test_Cheby_calore1_mm.xml
  test_Jacobi_nos1_hb.xml test_2_Jacobi_nos1_hb.xml nos1.rsa
  test_2_ILUT_nos1_hb.xml test_2_RILUK_nos1_hb.xml
  test_4_ILUT_nos1_hb.xml test_4_RILUK_nos1_hb.xml
  test_SGS_calore1_mm.xml calore1.mtx test_pseudoblockcg_small_sym_mm.xml
  test_Jacobi_small_sym_mm.xml small_sym.mtx test_gmres_small_sym_mm.xml
  test_Diagonal_calore1_mm.xml
  test_GS_tfqmr_small_sym_mm.xml test_ILUT_tfqmr_small_sym_mm.xml
  test_tfqmr_small_sym_mm.xml test_ILUT_tfqmr_calore1_mm.xml
  5w.mtx 5w.vec 5w_bel_tif_ILUT.xml
  5w_missing_diag.mtx 5w_missing_diag_ILUT.xml
  5w_bel_tif_RILUK_0.xml 5w_bel_tif_RILUK_1.xml
  test_RILUK_tfqmr_small_sym_mm.xml
  test_bordered_DIAG_small.xml small_sym_sing.mtx consist_rhs.vec small_null.vec


  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  EXEDEPS tif_belos
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME small_gmres_belos
  ARGS "--xml_file=test_gmres_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME small_pseudoblkcg_belos
  ARGS "--xml_file=test_pseudoblockcg_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Diagonal_belos
  ARGS "--xml_file=test_Diagonal_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Jac_sm_belos
  ARGS "--xml_file=test_Jacobi_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Jacobi_belos
  ARGS "--xml_file=test_Jacobi_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Jacobi_hb_belos
  ARGS "--xml_file=test_Jacobi_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Jacobi_hb_belos
  ARGS "--xml_file=test_2_Jacobi_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME ILUT_hb_belos
  ARGS "--xml_file=test_2_ILUT_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME ILUT_hb_belos
  ARGS "--xml_file=test_4_ILUT_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME RILUK_hb_belos
  ARGS "--xml_file=test_2_RILUK_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME RILUK_hb_belos
  ARGS "--xml_file=test_4_RILUK_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_TEST(
    tif_extprec_belos Ifpack2BelosCopyFiles
    NAME small_qd_tfqmr_belos
    ARGS "--xml_file=test_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_TEST(
    tif_extprec_belos Ifpack2BelosCopyFiles
    NAME small_qd_GS_tfqmr_belos
    ARGS "--xml_file=test_GS_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_TEST(
    tif_extprec_belos Ifpack2BelosCopyFiles
    NAME small_qd_ILUT_tfqmr_belos
    ARGS "--xml_file=test_ILUT_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_TEST(
    tif_extprec_belos Ifpack2BelosCopyFiles
    NAME qd_ILUT_tfqmr_belos
    ARGS "--xml_file=test_ILUT_tfqmr_calore1_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_TEST(
    tif_extprec_belos Ifpack2BelosCopyFiles
    NAME qd_ILUT_pbgmres_belos
    ARGS "--xml_file=5w_bel_tif_ILUT.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME ILUT_5w_2
  ARGS "--xml_file=5w_bel_tif_ILUT.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME ILUT_5w_no_diag
  ARGS "--xml_file=5w_missing_diag_ILUT.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

# mfh 31 Jan 2013: Disabled some of the RILUK tests for now, due to
# failures relating to sensitivity of the RILUK implementation to the
# Maps.  See my comments in the following (31 Jan 2013) commit:
#
# e486af8c46174efd16fb9e45792f8fee8e4d8859
#
# TRIBITS_ADD_TEST(
#   tif_belos Ifpack2BelosCopyFiles
#   NAME RILUK_5w_0
#   ARGS "--xml_file=5w_bel_tif_RILUK_0.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

# TRIBITS_ADD_TEST(
#   tif_belos Ifpack2BelosCopyFiles
#   NAME RILUK_5w_1
#   ARGS "--xml_file=5w_bel_tif_RILUK_1.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

# TRIBITS_ADD_TEST(
#   tif_belos Ifpack2BelosCopyFiles
#   NAME RILUK_small_belos
#   ARGS "--xml_file=test_RILUK_tfqmr_small_sym_mm.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Cheby_belos
  ARGS "--xml_file=test_Cheby_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME GS_belos
  ARGS "--xml_file=test_GS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME SGS_belos
  ARGS "--xml_file=test_SGS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME ILUT_belos
  ARGS "--xml_file=test_ILUT_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

# TRIBITS_ADD_TEST(
#   tif_belos Ifpack2BelosCopyFiles
#   NAME RILUK_belos
#   ARGS "--xml_file=test_RILUK_calore1_mm.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

TRIBITS_ADD_TEST(
  tif_belos Ifpack2BelosCopyFiles
  NAME Bordered_belos
  ARGS "--xml_file=test_bordered_DIAG_small.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

