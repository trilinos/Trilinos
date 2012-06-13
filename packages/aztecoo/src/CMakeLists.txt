
INCLUDE(TribitsLibraryMacros)

#
# A) Package-specific configuration options
#

SET(HAVE_BOOL ON)
SET(HAVE_STL ON)

IF (${PROJECT_NAME}_ENABLE_Fortran)
  SET(HAVE_FORTRAN_SUPPORT ON)
ENDIF()

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and directories)
#

#
# src
#

SET(HEADERS "")
SET(SOURCES "")

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS ${HEADERS}
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#
# C files
#

SET(HEADERS ${HEADERS}
  az_aztec_defs.h
  az_aztecf.h
  az_aztec.h
  az_f77func.h
  az_y12m_wrappers.h
  az_blas_wrappers.h
  az_lapack_wrappers.h
  )

SET(SOURCES ${SOURCES}
  az_comm.c 
  )

#
# C++ files
#

SET(HEADERS ${HEADERS}
  AztecOO.h
  AztecOO_string_maps.h
  AztecOO_Operator.h
  Epetra_MsrMatrix.h
  Aztec2Petra.h
  AZOO_iterate.h
  AZOO_printf.h
  AztecOO_Scaling.h
  AztecOO_StatusTest.h
  AztecOO_StatusType.h
  AztecOO_StatusTestCombo.h
  AztecOO_StatusTestMaxIters.h
  AztecOO_StatusTestResNorm.h
  AztecOO_ConditionNumber.h
  AztecOO_ConfigDefs.h
  )

SET(SOURCES ${SOURCES}
  AztecOO.cpp
  AztecOO_string_maps.cpp
  AztecOO_Operator.cpp
  Epetra_MsrMatrix.cpp
  Aztec2Petra.cpp
  AZOO_iterate.cpp
  AZOO_printf.cpp
  AztecOO_Scaling.cpp
  AztecOO_StatusTestCombo.cpp
  AztecOO_StatusTestMaxIters.cpp
  AztecOO_StatusTestResNorm.cpp
  AztecOO_ConditionNumber.cpp
  )

SET(HEADERS ${HEADERS}
  AztecOO_Version.h 
  )

SET(SOURCES ${SOURCES}
  az_check.c
  az_flop_cnt.c
  az_ilu_util.c
  az_interface.c
  az_matrix_util.c
  az_pad_utils.c
  az_tools.c
  az_util.c
  )

IF (TPL_ENABLE_MPI)

  SET(SOURCES ${SOURCES}
    md_timer_mpi.c
    md_wrap_mpi_c.c 
    )
ELSE()
  SET(HEADERS ${HEADERS}
    md_timer_generic.c
    md_wrap_scalar_c.c 
    )
ENDIF()

SET(SOURCES ${SOURCES}
  az_bilu.c
  az_cg.c
  az_cgs.c
  az_cgstab.c
  az_converge.c
  az_dd_overlap.c
  az_dgemv2.c
  az_dgemv3.c
  az_domain_decomp.c
  az_fix_pt.c
  az_gmres.c
  az_gmresr.c
  az_gsumd_puma.c
  az_icc.c
  az_ilut.c
  az_lu_y12.c
  az_poly.c
  az_printf.c
  az_qmrcgs.c
  az_rilu.c
  az_subdomain_solver.c
  az_cg_condnum.c
  az_gmres_condnum.c
  )

SET(SOURCES ${SOURCES}
  az_colperm.c
  az_matvec_mult.c
  az_old_matvec_mult.c
  az_precond.c
  az_scaling.c
  az_solve.c
  az_sort.c 
  )

IF (${PROJECT_NAME}_ENABLE_Fortran)

  SET(SOURCES ${SOURCES}
    az_reorder.f
    az_f_util.f
    )

ELSE()

  SET(SOURCES ${SOURCES}
  az_c_reorder.c
  az_c_util.c
  )

ENDIF()

#
# C) Define the targets for package's library(s)
#

TRIBITS_ADD_LIBRARY(
  aztecoo
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
