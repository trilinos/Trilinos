
INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

#SET(HAVE_BOOL ON)
#SET(HAVE_STL ON)

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

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

#
# support files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES ${SOURCES}
  Amesos_Status.cpp
  Amesos_Control.cpp
  Amesos_Support.cpp
  )

#
# DSCPACK files
#

IF (${PACKAGE_NAME}_ENABLE_DSCPACK)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Dscpack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Dscpack.cpp
    )
ENDIF()

#
# UMFPACK files
#

IF (${PACKAGE_NAME}_ENABLE_UMFPACK)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Umfpack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Umfpack.cpp
    )
ENDIF()

#
# MUMPS files
#

IF (${PACKAGE_NAME}_ENABLE_MUMPS)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Mumps.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Mumps.cpp
    )
ENDIF()

#
# PARDISO files
#

IF (${PACKAGE_NAME}_ENABLE_PARDISO)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Pardiso.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Pardiso.cpp
    )
ENDIF()

#
# TAUCS files
#

IF (${PACKAGE_NAME}_ENABLE_TAUCS)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Taucs.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Taucs.cpp
    )
ENDIF()

#
# PARAKLETE KLU files
#

SET(PARAKLETE_KLU_COMMON_H
  SuiteSparse/AMD/Include/amesos_amd.h
  SuiteSparse/AMD/Include/amesos_amd_internal.h
  SuiteSparse/COLAMD/Include/amesos_colamd.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod_check.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod_core.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod_internal.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod_io64.h
  SuiteSparse/CHOLMOD/Include/amesos_cholmod_partition.h
  SuiteSparse/UFconfig/amesos_UFconfig.h
  )

SET(PARAKLETE_KLU_COMMON
  SuiteSparse/AMD/Source/amesos_amd_1.c
  SuiteSparse/AMD/Source/amesos_amd_2.c
  SuiteSparse/AMD/Source/amesos_amd_aat.c
  SuiteSparse/AMD/Source/amesos_amd_control.c
  SuiteSparse/AMD/Source/amesos_amd_defaults.c
  SuiteSparse/AMD/Source/amesos_amd_dump.c
  SuiteSparse/AMD/Source/amesos_amd_global.c
  SuiteSparse/AMD/Source/amesos_amd_info.c
  SuiteSparse/AMD/Source/amesos_amd_order.c
  SuiteSparse/AMD/Source/amesos_amd_postorder.c
  SuiteSparse/AMD/Source/amesos_amd_post_tree.c
  SuiteSparse/AMD/Source/amesos_amd_preprocess.c
  SuiteSparse/AMD/Source/amesos_amd_valid.c
  SuiteSparse/COLAMD/Source/amesos_colamd.c
  SuiteSparse/COLAMD/Source/amesos_colamd_global.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_common.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_complex.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_error.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_memory.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_sparse.c
  SuiteSparse/CHOLMOD/Core/amesos_cholmod_transpose.c
  SuiteSparse/CAMD/Source/amesos_camd_2.c
  SuiteSparse/CAMD/Source/amesos_camd_postorder.c
  )

#
# KLU files
#

IF (${PACKAGE_NAME}_ENABLE_KLU)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/AMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/COLAMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/UFconfig)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/KLU/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/BTF/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CAMD/Include)

  SET(HEADERS ${HEADERS}
    ${PARAKLETE_KLU_COMMON_H}
    SuiteSparse/KLU/Include/amesos_klu_decl.h
    SuiteSparse/KLU/Include/amesos_klu_internal.h
    SuiteSparse/KLU/Include/amesos_klu_version.h
    SuiteSparse/BTF/Include/amesos_btf_decl.h
    SuiteSparse/BTF/Include/amesos_btf_internal.h
    Amesos_Klu.h
    )

  SET(SOURCES ${SOURCES}
    ${PARAKLETE_KLU_COMMON}
    Amesos_Klu.cpp
    SuiteSparse/KLU/Source/amesos_klu_analyze.c
    SuiteSparse/KLU/Source/amesos_klu_analyze_given.c
    SuiteSparse/KLU/Source/amesos_klu_impl.c
    SuiteSparse/KLU/Source/amesos_klu_defaults.c
    SuiteSparse/KLU/Source/amesos_klu_diagnostics.c
    SuiteSparse/KLU/Source/amesos_klu_dump.c
    SuiteSparse/KLU/Source/amesos_klu_extract.c
    SuiteSparse/KLU/Source/amesos_klu_factor.c
    SuiteSparse/KLU/Source/amesos_klu_free_numeric.c
    SuiteSparse/KLU/Source/amesos_klu_free_symbolic.c
    SuiteSparse/KLU/Source/amesos_klu_kernel.c
    SuiteSparse/KLU/Source/amesos_klu_memory.c
    SuiteSparse/KLU/Source/amesos_klu_refactor.c
    SuiteSparse/KLU/Source/amesos_klu_scale.c
    SuiteSparse/KLU/Source/amesos_klu_solve.c
    SuiteSparse/KLU/Source/amesos_klu_sort.c
    SuiteSparse/KLU/Source/amesos_klu_tsolve.c
    SuiteSparse/BTF/Source/amesos_btf_order.c
    SuiteSparse/BTF/Source/amesos_btf_maxtrans.c
    SuiteSparse/BTF/Source/amesos_btf_strongcomp.c
    )
ENDIF()

#
# PARAKLETE files
#

IF (${PACKAGE_NAME}_ENABLE_PARAKLETE)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/AMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/COLAMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CHOLMOD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/UFconfig)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/PARAKLETE/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CAMD/Include)

  SET(HEADERS ${HEADERS}
    ${PARAKLETE_KLU_COMMON_H}
    Amesos_Paraklete.h
    SuiteSparse/PARAKLETE/Include/amesos_paraklete_decl.h
    SuiteSparse/CHOLMOD/Include/amesos_cholmod_complexity.h
    SuiteSparse/CCOLAMD/Include/amesos_ccolamd.h
    SuiteSparse/CHOLMOD/Include/amesos_cholmod_blas.h
    SuiteSparse/CHOLMOD/Include/amesos_cholmod_cholesky.h
    SuiteSparse/CHOLMOD/Include/amesos_cholmod_config.h
    SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_dense.c
    SuiteSparse/CHOLMOD/Include/amesos_cholmod_template.h
    SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_change_factor.c
    SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_transpose.c
    SuiteSparse/CHOLMOD/Core/amesos_t_cholmod_triplet.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_solve.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_lsolve.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_t_cholmod_ltsolve.c
    SuiteSparse/CAMD/Include/amesos_camd.h
    SuiteSparse/CAMD/Include/amesos_camd_internal.h
    )

  SET(SOURCES ${SOURCES}
    ${PARAKLETE_KLU_COMMON}
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_analyze.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_btf.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_error.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_factorize.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_finish.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_kernel.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_lsolve_node.c 
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_node.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_reanalyze.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_solve.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_start.c
    SuiteSparse/PARAKLETE/Source/amesos_paraklete_usolve_node.c 
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_common.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_dense.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_factor.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_change_factor.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_memory.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_complex.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_transpose.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_sparse.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_band.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_copy.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_triplet.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_error.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_aat.c
    SuiteSparse/CHOLMOD/Core/amesos_cholmod_l_add.c
    SuiteSparse/CHOLMOD/Check/amesos_cholmod_check.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_postorder.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_amd.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_colamd.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_etree.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_rowcolcounts.c
    SuiteSparse/CHOLMOD/Cholesky/amesos_cholmod_l_analyze.c
    SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_ccolamd.c
    SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_camd.c
    SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_csymamd.c
    SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_metis.c
    SuiteSparse/CHOLMOD/Partition/amesos_cholmod_l_nesdis.c
    SuiteSparse/AMD/Source/amesos_amd_l1.c
    SuiteSparse/AMD/Source/amesos_amd_l2.c
    SuiteSparse/AMD/Source/amesos_amd_l_aat.c
    SuiteSparse/AMD/Source/amesos_amd_l_order.c
    SuiteSparse/AMD/Source/amesos_amd_l_postorder.c
    SuiteSparse/AMD/Source/amesos_amd_l_post_tree.c
    SuiteSparse/AMD/Source/amesos_amd_l_preprocess.c
    SuiteSparse/AMD/Source/amesos_amd_l_valid.c
    SuiteSparse/CAMD/Source/amesos_camd_global.c
    SuiteSparse/CAMD/Source/amesos_camd_l2.c
    SuiteSparse/CAMD/Source/amesos_camd_l_postorder.c
    SuiteSparse/COLAMD/Source/amesos_colamd_l.c
    SuiteSparse/CCOLAMD/Source/amesos_ccolamd_l.c
    SuiteSparse/CCOLAMD/Source/amesos_ccolamd_global.c
    SuiteSparse/KLU/Source/amesos_klu_l.c
    SuiteSparse/KLU/Source/amesos_klu_l_analyze.c
    SuiteSparse/KLU/Source/amesos_klu_l_analyze_given.c
    SuiteSparse/KLU/Source/amesos_klu_l_defaults.c
    SuiteSparse/KLU/Source/amesos_klu_l_free_numeric.c
    SuiteSparse/KLU/Source/amesos_klu_l_free_symbolic.c
    SuiteSparse/KLU/Source/amesos_klu_l_factor.c
    SuiteSparse/KLU/Source/amesos_klu_l_kernel.c
    SuiteSparse/KLU/Source/amesos_klu_l_memory.c
    SuiteSparse/KLU/Source/amesos_klu_l_scale.c
    SuiteSparse/KLU/Source/amesos_klu_l_solve.c
    SuiteSparse/BTF/Source/amesos_btf_l_order.c
    SuiteSparse/BTF/Source/amesos_btf_l_maxtrans.c
    SuiteSparse/BTF/Source/amesos_btf_l_strongcomp.c
    Amesos_Paraklete.cpp
    )
ENDIF()

#
# MC64 files
#

IF (${PACKAGE_NAME}_ENABLE_MC64)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_MC64.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_MC64.cpp
    )
ENDIF()

#
# SCALAPACK files
#

IF (${PACKAGE_NAME}_ENABLE_SCALAPACK)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Scalapack.h
    Amesos_SCALAPACK_wrappers.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Scalapack.cpp
    )
ENDIF()

#
# FACTORY files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Amesos.h
  )

SET(SOURCES ${SOURCES}
  Amesos.cpp
  )

#
# SUPERLU files
#

IF (${PACKAGE_NAME}_ENABLE_SUPERLU)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Superlu.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Superlu.cpp
    )
ENDIF()

#
# SUPERLUDIST files
#

IF (${PACKAGE_NAME}_ENABLE_SuperLUDist)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Superludist.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Superludist.cpp
    )
ENDIF()

#
# COMP files
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Amesos_BaseSolver.h
  Amesos_ComponentBaseSolver.h
  Amesos_ConfigDefs.h
  CrsMatrixTranspose.h
  Amesos_NoCopiable.h
  Amesos_Utils.h
  Amesos_Time.h
  Amesos_Control.h
  Amesos_Support.h
  Amesos_Status.h
  Amesos_Scaling.h
  Amesos_Reordering.h
  )

SET(SOURCES ${SOURCES}
  CrsMatrixTranspose.cpp  
  )

IF (${PACKAGE_NAME}_ENABLE_LAPACK)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Lapack.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Lapack.cpp
    )
ENDIF()



##############################

#
# src
#

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS ${HEADERS}
  Amesos_BaseSolver.h
  Amesos_ComponentBaseSolver.h
  Amesos_TestRowMatrix.h
  )

SET(SOURCES ${SOURCES}

  )

IF (TPL_ENABLE_MPI)

  SET(SOURCES ${SOURCES}

    )

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  amesos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
