
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

ADD_SUBDIRECTORY(AdditiveSchwarzRILUK)

IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE(
    tif_belos
    CATEGORIES BASIC PERFORMANCE
    SOURCES belos_solve.cpp
    )
ENDIF()

IF(Ifpack2_ENABLE_QD)
  TRIBITS_ADD_EXECUTABLE(
    tif_extprec_belos
    SOURCES belos_extprec_solve.cpp
  )
ENDIF()

IF (Tpetra_INST_LONG_DOUBLE)
  TRIBITS_ADD_EXECUTABLE(
    tif_belos_longdouble
    SOURCES belos_solve_longdouble.cpp
    )
ENDIF()

IF (Tpetra_ENABLE_quadmath)
  TRIBITS_ADD_EXECUTABLE(
    tif_belos_float128
    SOURCES belos_solve_float128.cpp
    )
ENDIF()

IF (Tpetra_INST_DOUBLE)
  TRIBITS_ADD_EXECUTABLE(
    tif_tpetra_native
    SOURCES tpetra_native.cpp
  )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(Ifpack2BelosCopyFiles
  DEST_FILES
  test_database_schwarz.xml
  test_ILUT_calore1_mm.xml
  test_Jacobi_calore1_mm.xml
  test_Jacobi_calore1_mm_constGraph.xml
  test_BlockRelaxationAmesos2_calore1_mm.xml
  test_BlockRelaxationZoltan2_calore1_mm.xml
  test_GS_calore1_mm.xml
  test_MTGS_calore1_mm.xml
  test_RILUK_calore1_mm.xml
  test_Cheby_calore1_mm.xml
  test_Cheby_calore1_nospectralradius_mm.xml
  test_FILU_calore1_mm.xml
  test_FILU_calore1_mm_sptrsv.xml
  test_FIC_calore1_mm.xml
  test_FIC_calore1_mm_sptrsv.xml
  test_FILDL_calore1_mm.xml
  test_FILDL_calore1_mm_sptrsv.xml
  test_Jacobi_nos1_hb.xml
  test_2_Jacobi_nos1_hb.xml
  test_Jacobi_gcrodr_nos1_hb.xml
  nos1.rsa
  test_Jacobi_bcsstk14_hb.xml
  test_Jacobi_minres_bcsstk14_hb.xml
  test_Jacobi_pseudoblockcg_bcsstk14_hb.xml
  bcsstk14.hb
  test_2_ILUT_nos1_hb.xml
  test_2_RILUK_nos1_hb.xml
  test_RBILUK_nos1_hb_block.xml
  test_2_RILUK_HTS_nos1_hb.xml
  test_2_RILUK_2streams_nos1_hb.xml
  test_2_RILUK_4streams_nos1_hb.xml
  test_2_RILUK_2streams_rcm_nos1_hb.xml
  test_2_RILUK_4streams_rcm_nos1_hb.xml
  test_4_ILUT_nos1_hb.xml
  test_4_RILUK_nos1_hb.xml
  test_4_RILUK_HTS_nos1_hb.xml
  test_4_RILUK_2streams_nos1_hb.xml
  test_4_RILUK_4streams_nos1_hb.xml
  test_4_RILUK_2streams_rcm_nos1_hb.xml
  test_4_RILUK_4streams_rcm_nos1_hb.xml
  test_SGS_calore1_mm.xml
  test_MTSGS_calore1_mm.xml
  calore1.mtx
  calore1_rhs.mtx
  test_pseudoblockcg_small_sym_mm.xml
  test_minres_small_sym_mm.xml
  test_Jacobi_small_sym_mm.xml
  small_sym.mtx
  test_gmres_small_sym_mm.xml
  test_Diagonal_tfqmr_calore1_mm.xml
  test_Diagonal_bicgstab_calore1_mm.xml
  test_Diagonal_gcrodr_calore1_mm.xml
  test_GS_tfqmr_small_sym_mm.xml
  test_ILUT_tfqmr_small_sym_mm.xml
  test_tfqmr_small_sym_mm.xml
  test_FILU_tfqmr_small_sym_mm.xml
  test_FILU_tfqmr_small_sym_mm_sptrsv.xml
  test_FILU_tfqmr_small_sym_mm_schwarz.xml
  test_FIC_tfqmr_small_sym_mm.xml
  test_FIC_tfqmr_small_sym_mm_sptrsv.xml
  test_FIC_tfqmr_small_sym_mm_schwarz.xml
  test_FILDL_tfqmr_small_sym_mm.xml
  test_FILDL_tfqmr_small_sym_mm_sptrsv.xml
  test_FILDL_tfqmr_small_sym_mm_schwarz.xml
  test_ILUT_tfqmr_calore1_mm.xml
  5w.mtx
  5w.vec
  6w.mtx
  6w.vec
  5w_bel_tif_ILUT.xml
  5w_bel_tif_FILU.xml
  6w_bel_tif_FILU_blockcrs.xml
  5w_bel_tif_FILU_sptrsv.xml
  5w_bel_tif_FILU_sptrsv_host.xml
  5w_bel_tif_FILU_schwarz.xml
  5w_bel_tif_FIC.xml
  5w_bel_tif_FIC_sptrsv.xml
  5w_bel_tif_FIC_schwarz.xml
  5w_bel_tif_FILDL.xml
  5w_bel_tif_FILDL_sptrsv.xml
  5w_bel_tif_FILDL_schwarz.xml
  5w_missing_diag.mtx
  5w_missing_diag_ILUT.xml
  5w_bel_tif_RILUK_0.xml
  5w_bel_tif_RILUK_1.xml
  test_RILUK_tfqmr_small_sym_mm.xml
  test_bordered_DIAG_small.xml
  test_minres_bordered_DIAG_small.xml
  test_FILU_small_sym_sing.xml
  test_FILU_small_sym_sing_sptrsv.xml
  test_FIC_small_sym_sing.xml
  test_FIC_small_sym_sing_sptrsv.xml
  test_FILDL_small_sym_sing.xml
  test_FILDL_small_sym_sing_sptrsv.xml
  small_sym_sing.mtx consist_rhs.vec small_null.vec
  mueluMatrix.mtx
  mueluMatrix_rhs.mtx
  test_4_GS_muelu.xml
  test_4_MTGS_muelu.xml
  test_4_SGS_muelu.xml
  test_4_MTSGS_muelu.xml
  schwarz.xml
  # test_SGS_wathen_mm.xml
  # test_MTSGS_wathen_mm.xml

  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  #EXEDEPS tif_belos
)

IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_TEST(
  tif_belos
  NAME small_gmres_belos
  ARGS "--xml_file=test_gmres_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
ENDIF()

IF (Tpetra_INST_LONG_DOUBLE)
TRIBITS_ADD_TEST(
  tif_belos_longdouble
  NAME gmres_belos_longdouble
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
ENDIF()

IF (Tpetra_ENABLE_quadmath)
TRIBITS_ADD_TEST(
  tif_belos_float128
  NAME gmres_belos_float128
  ARGS ""
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
ENDIF()


IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_TEST(
  tif_belos
  NAME small_pseudoblkcg_belos
  ARGS "--xml_file=test_pseudoblockcg_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME small_minres_belos
  ARGS "--xml_file=test_minres_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Diagonal_tfqmr_belos
  ARGS "--xml_file=test_Diagonal_tfqmr_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Diagonal_bicgstab_belos
  ARGS "--xml_file=test_Diagonal_bicgstab_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Diagonal_gcrodr_belos
  ARGS "--xml_file=test_Diagonal_gcrodr_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jac_sm_belos
  ARGS "--xml_file=test_Jacobi_small_sym_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_belos
  ARGS "--xml_file=test_Jacobi_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME DatabaseSchwarz_belos
  ARGS "--xml_file=test_database_schwarz.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_belos_constGraph
  ARGS "--xml_file=test_Jacobi_calore1_mm_constGraph.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1-4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_hb_belos
  ARGS "--xml_file=test_Jacobi_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_gcrodr_hb_belos
  ARGS "--xml_file=test_Jacobi_gcrodr_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_hb_belos
  ARGS "--xml_file=test_2_Jacobi_nos1_hb.xml"
  COMM serial mpi NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

# HKT: Commenting out test temporarily because it is causing random failures
#TRIBITS_ADD_TEST(
#  tif_belos
#  NAME Jacobi_bcsstk14_hb_belos
#  ARGS "--xml_file=test_Jacobi_bcsstk14_hb.xml"
#  COMM serial mpi
#  NUM_MPI_PROCS 1
#  STANDARD_PASS_OUTPUT
#)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_bcsstk14_minres_hb_belos
  ARGS "--xml_file=test_Jacobi_minres_bcsstk14_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Jacobi_bcsstk14_pbcg_belos
  ARGS "--xml_file=test_Jacobi_pseudoblockcg_bcsstk14_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME ILUT_hb_belos
  ARGS "--xml_file=test_2_ILUT_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME ILUT_hb_belos
  ARGS "--xml_file=test_4_ILUT_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_hb_belos
  ARGS "--xml_file=test_2_RILUK_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RBILUK_hb_belos_block_serial
  ARGS "--xml_file=test_RBILUK_nos1_hb_block.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_hb_belos
  ARGS "--xml_file=test_4_RILUK_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core AND ${PACKAGE_NAME}_ENABLE_Xpetra)
  TRIBITS_ADD_TEST(
    tif_belos
    NAME Schwarz_RILUK_hb_belos
    ARGS "--xml_file=schwarz.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
  )
ENDIF()




TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_2streams_hb_belos
  ARGS "--xml_file=test_2_RILUK_2streams_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_4streams_hb_belos
  ARGS "--xml_file=test_2_RILUK_4streams_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_2streams_hb_belos
  ARGS "--xml_file=test_4_RILUK_2streams_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_4streams_hb_belos
  ARGS "--xml_file=test_4_RILUK_4streams_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)
TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_2streams_rcm_hb_belos
  ARGS "--xml_file=test_2_RILUK_2streams_rcm_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_4streams_rcm_hb_belos
  ARGS "--xml_file=test_2_RILUK_4streams_rcm_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_2streams_rcm_hb_belos
  ARGS "--xml_file=test_4_RILUK_2streams_rcm_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME RILUK_4streams_rcm_hb_belos
  ARGS "--xml_file=test_4_RILUK_4streams_rcm_nos1_hb.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)
ENDIF()


IF(Tpetra_INST_DOUBLE AND Ifpack2_ENABLE_ShyLU_NodeHTS)
  TRIBITS_ADD_TEST(
    tif_belos
    NAME RILUK_HTS_hb_belos
    ARGS "--xml_file=test_2_RILUK_HTS_nos1_hb.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    STANDARD_PASS_OUTPUT
    )

  TRIBITS_ADD_TEST(
    tif_belos
    NAME RILUK_HTS_hb_belos
    ARGS "--xml_file=test_4_RILUK_HTS_nos1_hb.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
    )
ENDIF()

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

IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_TEST(
  tif_belos
  NAME ILUT_5w_2
  ARGS "--xml_file=5w_bel_tif_ILUT.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME ILUT_5w_no_diag
  ARGS "--xml_file=5w_missing_diag_ILUT.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
ENDIF()

IF(Tpetra_INST_DOUBLE AND Ifpack2_ENABLE_ShyLU_NodeFastILU)
  #calore1
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_calore1
    ARGS "--xml_file=test_FILU_calore1_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )

  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_calore1_sptrsv
    ARGS "--xml_file=test_FILU_calore1_mm_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_calore1
    ARGS "--xml_file=test_FILDL_calore1_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_calore1_sptrsv
    ARGS "--xml_file=test_FILDL_calore1_mm_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  #5w
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_5w
    ARGS "--xml_file=5w_bel_tif_FILU.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_6w_blockcrs
    ARGS "--xml_file=6w_bel_tif_FILU_blockcrs.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_5w_sptrsv
    ARGS "--xml_file=5w_bel_tif_FILU_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_5w_sptrsv_host
    ARGS "--xml_file=5w_bel_tif_FILU_sptrsv_host.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FILU_5w_schwarz
      ARGS "--xml_file=5w_bel_tif_FILU_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_5w
    ARGS "--xml_file=5w_bel_tif_FIC.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_5w_sptrsv
    ARGS "--xml_file=5w_bel_tif_FIC_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FIC_5w_schwarz
      ARGS "--xml_file=5w_bel_tif_FIC_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_5w
    ARGS "--xml_file=5w_bel_tif_FILDL.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_5w_sptrsv
    ARGS "--xml_file=5w_bel_tif_FILDL_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FILDL_5w_schwarz
      ARGS "--xml_file=5w_bel_tif_FILDL_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  #small_sym
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_small_sym
    ARGS "--xml_file=test_FILU_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_small_sym_sptrsv
    ARGS "--xml_file=test_FILU_tfqmr_small_sym_mm_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FILU_small_sym_schwarz
      ARGS "--xml_file=test_FILU_tfqmr_small_sym_mm_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_small_sym
    ARGS "--xml_file=test_FIC_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_small_sym_sptrsv
    ARGS "--xml_file=test_FIC_tfqmr_small_sym_mm_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FIC_small_sym_schwarz
      ARGS "--xml_file=test_FIC_tfqmr_small_sym_mm_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_small_sym
    ARGS "--xml_file=test_FILDL_tfqmr_small_sym_mm.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_small_sym_sptrsv
    ARGS "--xml_file=test_FILDL_tfqmr_small_sym_mm_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    TRIBITS_ADD_TEST(
      tif_belos
      NAME FILDL_small_sym_schwarz
      ARGS "--xml_file=test_FILDL_tfqmr_small_sym_mm_schwarz.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      STANDARD_PASS_OUTPUT
    )
  ENDIF()
  #small_sym_sing
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_small_sym_sing
    ARGS "--xml_file=test_FILU_small_sym_sing.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILU_small_sym_sing_sptrsv
    ARGS "--xml_file=test_FILU_small_sym_sing_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_small_sym_sing
    ARGS "--xml_file=test_FIC_small_sym_sing.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FIC_small_sym_sing_sptrsv
    ARGS "--xml_file=test_FIC_small_sym_sing_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_small_sym_sing
    ARGS "--xml_file=test_FILDL_small_sym_sing.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
  TRIBITS_ADD_TEST(
    tif_belos
    NAME FILDL_small_sym_sing_sptrsv
    ARGS "--xml_file=test_FILDL_small_sym_sing_sptrsv.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
  )
ENDIF()

# mfh 31 Jan 2013: Disabled some of the RILUK tests for now, due to
# failures relating to sensitivity of the RILUK implementation to the
# Maps.  See my comments in the following (31 Jan 2013) commit:
#
# e486af8c46174efd16fb9e45792f8fee8e4d8859
#
# TRIBITS_ADD_TEST(
#   tif_belos
#   NAME RILUK_5w_0
#   ARGS "--xml_file=5w_bel_tif_RILUK_0.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

# TRIBITS_ADD_TEST(
#   tif_belos
#   NAME RILUK_5w_1
#   ARGS "--xml_file=5w_bel_tif_RILUK_1.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

# TRIBITS_ADD_TEST(
#   tif_belos
#   NAME RILUK_small_belos
#   ARGS "--xml_file=test_RILUK_tfqmr_small_sym_mm.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )

IF (Tpetra_INST_DOUBLE)
TRIBITS_ADD_TEST(
  tif_belos
  NAME Cheby_belos
  ARGS "--xml_file=test_Cheby_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

# GH 04 Nov 2022: Test the non-spectral
# radius route of the power method for robustness
TRIBITS_ADD_TEST(
  tif_belos
  NAME Cheby_belos_nospectralradius
  ARGS "--xml_file=test_Cheby_calore1_nospectralradius_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME GS_belos
  ARGS "--xml_file=test_GS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME SGS_belos
  ARGS "--xml_file=test_SGS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME MTGS_belos
  ARGS "--xml_file=test_MTGS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME MTSGS_belos
  ARGS "--xml_file=test_MTSGS_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME GS_belos_muelu
  ARGS "--xml_file=test_4_GS_muelu.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME MTGS_belos_muelu
  ARGS "--xml_file=test_4_MTGS_muelu.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME SGS_belos_muelu
  ARGS "--xml_file=test_4_SGS_muelu.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )

TRIBITS_ADD_TEST(
  tif_belos
  NAME MTSGS_belos_muelu
  ARGS "--xml_file=test_4_MTSGS_muelu.xml"
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
  )


TRIBITS_ADD_TEST(
  tif_belos
  NAME ILUT_belos
  ARGS "--xml_file=test_ILUT_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

# TRIBITS_ADD_TEST(
#   tif_belos
#   NAME RILUK_belos
#   ARGS "--xml_file=test_RILUK_calore1_mm.xml"
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   STANDARD_PASS_OUTPUT
# )


TRIBITS_ADD_TEST(
  tif_belos
  NAME Bordered_blockcg_belos
  ARGS "--xml_file=test_bordered_DIAG_small.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_belos
  NAME Bordered_minres_belos
  ARGS "--xml_file=test_minres_bordered_DIAG_small.xml"
  COMM serial mpi
  NUM_MPI_PROCS 1
  STANDARD_PASS_OUTPUT
)
ENDIF()

TRIBITS_ADD_TEST(
  tif_tpetra_native
  NAME Tpetra_Native_GMRES_SINGLE_REDUCE
  ARGS "--matrixFilename=calore1.mtx --convergenceTolerances=1e-8 --solverTypes=\"TPETRA GMRES SINGLE REDUCE\""
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_tpetra_native
  NAME Tpetra_Native_GMRES_PIPELINE
  ARGS "--matrixFilename=calore1.mtx --convergenceTolerances=1e-8 --solverTypes=\"TPETRA GMRES PIPELINE\""
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

TRIBITS_ADD_TEST(
  tif_tpetra_native
  NAME Tpetra_Native_GMRES_S_STEP
  ARGS "--matrixFilename=calore1.mtx --convergenceTolerances=1e-8 --solverTypes=\"TPETRA GMRES S-STEP\""
  COMM serial mpi
  NUM_MPI_PROCS 4
  STANDARD_PASS_OUTPUT
)

IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
TRIBITS_ADD_TEST(
  tif_belos
  NAME BlockRelaxation_Zoltan2_belos
  ARGS "--xml_file=test_BlockRelaxationZoltan2_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Amesos2)
TRIBITS_ADD_TEST(
  tif_belos
  NAME BlockRelaxation_Amesos2_belos
  ARGS "--xml_file=test_BlockRelaxationAmesos2_calore1_mm.xml"
  COMM serial mpi
  NUM_MPI_PROCS 2
  STANDARD_PASS_OUTPUT
)
ENDIF()
