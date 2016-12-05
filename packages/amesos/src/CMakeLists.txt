

#
# A) Package-specific configuration options
#

#SET(HAVE_BOOL ON)
#SET(HAVE_STL ON)

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
   
   IF(TPL_ENABLE_MPI)

     INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

     SET(HEADERS ${HEADERS}
     Amesos_Mumps.h
     )

     SET(SOURCES ${SOURCES}
     Amesos_Mumps.cpp
     )
  ELSE()

     MESSAGE(WARNING "\n\n***Amesos with MUMPS requires MPI, MUMPS will not be available***\n\n")
   
   ENDIF()     

ENDIF()

#
# PARDISO files
#

IF (${PACKAGE_NAME}_ENABLE_PARDISO OR ${PACKAGE_NAME}_ENABLE_PARDISO_MKL)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_Pardiso.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_Pardiso.cpp
    )
ENDIF()

#
# CSparse files
#

IF (${PACKAGE_NAME}_ENABLE_CSparse OR ${PACKAGE_NAME}_ENABLE_CSparse)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

  SET(HEADERS ${HEADERS}
    Amesos_CSparse.h
    )

  SET(SOURCES ${SOURCES}
    Amesos_CSparse.cpp
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

SET(TRILINOS_KLU_COMMON_H
  SuiteSparse/AMD/Include/trilinos_amd.h
  SuiteSparse/AMD/Include/trilinos_amd_internal.h
  SuiteSparse/COLAMD/Include/trilinos_colamd.h
  SuiteSparse/UFconfig/amesos_UFconfig.h
  )

SET(TRILINOS_KLU_COMMON
  SuiteSparse/AMD/Source/trilinos_amd_1.c
  SuiteSparse/AMD/Source/trilinos_amd_2.c
  SuiteSparse/AMD/Source/trilinos_amd_aat.c
  SuiteSparse/AMD/Source/trilinos_amd_control.c
  SuiteSparse/AMD/Source/trilinos_amd_defaults.c
  SuiteSparse/AMD/Source/trilinos_amd_dump.c
  SuiteSparse/AMD/Source/trilinos_amd_global.c
  SuiteSparse/AMD/Source/trilinos_amd_info.c
  SuiteSparse/AMD/Source/trilinos_amd_order.c
  SuiteSparse/AMD/Source/trilinos_amd_postorder.c
  SuiteSparse/AMD/Source/trilinos_amd_post_tree.c
  SuiteSparse/AMD/Source/trilinos_amd_preprocess.c
  SuiteSparse/AMD/Source/trilinos_amd_valid.c
  SuiteSparse/AMD/Source/trilinos_amd_l1.c
  SuiteSparse/AMD/Source/trilinos_amd_l2.c
  SuiteSparse/AMD/Source/trilinos_amd_l_aat.c
  SuiteSparse/AMD/Source/trilinos_amd_l_order.c
  SuiteSparse/AMD/Source/trilinos_amd_l_postorder.c
  SuiteSparse/AMD/Source/trilinos_amd_l_post_tree.c
  SuiteSparse/AMD/Source/trilinos_amd_l_preprocess.c
  SuiteSparse/AMD/Source/trilinos_amd_l_valid.c
  SuiteSparse/COLAMD/Source/trilinos_colamd.c
  SuiteSparse/COLAMD/Source/trilinos_colamd_l.c
  SuiteSparse/COLAMD/Source/trilinos_colamd_global.c
  SuiteSparse/CAMD/Source/trilinos_camd_1.c
  SuiteSparse/CAMD/Source/trilinos_camd_defaults.c
  SuiteSparse/CAMD/Source/trilinos_camd_l2.c
  SuiteSparse/CAMD/Source/trilinos_camd_preprocess.c
  SuiteSparse/CAMD/Source/trilinos_camd_2.c
  SuiteSparse/CAMD/Source/trilinos_camd_dump.c
  SuiteSparse/CAMD/Source/trilinos_camd_l_postorder.c
  SuiteSparse/CAMD/Source/trilinos_camd_valid.c
  SuiteSparse/CAMD/Source/trilinos_camd_aat.c
  SuiteSparse/CAMD/Source/trilinos_camd_global.c
  SuiteSparse/CAMD/Source/trilinos_camd_order.c
  SuiteSparse/CAMD/Source/trilinos_camd_control.c
  SuiteSparse/CAMD/Source/trilinos_camd_info.c
  SuiteSparse/CAMD/Source/trilinos_camd_postorder.c
  )


  IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
      IF(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.2)
        SET_SOURCE_FILES_PROPERTIES(${TRILINOS_KLU_COMMON},
                PROPERTIES
                COMPILE_FLAGS "-Wno-all -Wno-pedantic")
      ENDIF()
  ELSEIF (WIN32 AND NOT CYGWIN)
       # No flags
  ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
      #No flags
  ELSE()
       SET_SOURCE_FILES_PROPERTIES(${TRILINOS_KLU_COMMON},
                PROPERTIES
                COMPILE_FLAGS "-Wno-all")
  ENDIF()
         

#
# KLU files
#

IF (${PACKAGE_NAME}_ENABLE_KLU)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/AMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/COLAMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/UFconfig)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/KLU/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/BTF/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CAMD/Include)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/SuiteSparse/CCOLAMD/Include)

  SET(HEADERS ${HEADERS}
    ${TRILINOS_KLU_COMMON_H}
    SuiteSparse/KLU/Include/amesos_klu_decl.h
    SuiteSparse/KLU/Include/amesos_klu_internal.h
    SuiteSparse/KLU/Include/amesos_klu_version.h
    SuiteSparse/BTF/Include/trilinos_btf_decl.h
    SuiteSparse/BTF/Include/trilinos_btf_internal.h
    SuiteSparse/CCOLAMD/Include/trilinos_ccolamd.h
    SuiteSparse/CAMD/Include/trilinos_camd.h
    SuiteSparse/CAMD/Include/trilinos_camd_internal.h
    Amesos_Klu.h
    )

   SET(KLU_FILES
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
    SuiteSparse/BTF/Source/trilinos_btf_order.c
    SuiteSparse/BTF/Source/trilinos_btf_maxtrans.c
    SuiteSparse/BTF/Source/trilinos_btf_strongcomp.c
    SuiteSparse/BTF/Source/trilinos_btf_l_order.c
    SuiteSparse/BTF/Source/trilinos_btf_l_maxtrans.c
    SuiteSparse/BTF/Source/trilinos_btf_l_strongcomp.c
    SuiteSparse/CCOLAMD/Source/trilinos_ccolamd.c
    SuiteSparse/CCOLAMD/Source/trilinos_ccolamd_l.c
    SuiteSparse/CCOLAMD/Source/trilinos_ccolamd_global.c
    )
    

    IF(CMAKE_COMPILER_IS_GNUCC OR CMAKE_COMPILER_IS_GNUCXX)
        IF(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 4.2)
        SET_SOURCE_FILES_PROPERTIES(${KLU_FILES},
                PROPERTIES
                COMPILE_FLAGS "-Wno-all -Wno-pedantic")
        ENDIF()
    ELSEIF (WIN32 AND NOT CYGWIN)
       # No flags
    ELSEIF ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "PGI")
      #No flags
    ELSE()
         SET_SOURCE_FILES_PROPERTIES(${KLU_FILES},
                PROPERTIES
                COMPILE_FLAGS "-Wno-all")
    ENDIF()



    SET(SOURCES ${SOURCES}
    ${TRILINOS_KLU_COMMON}
    Amesos_Klu.cpp
    ${KLU_FILES}
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

IF (${PACKAGE_NAME}_ENABLE_SuperLU)
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

TRIBITS_ADD_LIBRARY(
  amesos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
