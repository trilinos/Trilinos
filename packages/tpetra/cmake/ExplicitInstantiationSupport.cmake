include(Join)

# Tpetra ETI type fields
SET(Tpetra_ETI_FIELDS "SIN|SOUT|S|LO|GO|N|CS|DS")

# Exclude all of the types that CUDA/Thrust doesn't support
ASSERT_DEFINED(Tpetra_ENABLE_Thrust)
IF(Tpetra_ENABLE_Thrust)
  # no dd_real/qd_real support for CUDA, nor int/complex even via Cusp :( 
  SET(CUDA_UNSUPPORTED_SCALARS "long|dd_real|qd_real|std::complex<double>|std::complex<float>")
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET     "S=${CUDA_UNSUPPORTED_SCALARS}"                 "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET     "SIN=.*" "SOUT=${CUDA_UNSUPPORTED_SCALARS}|int" "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET     "SIN=${CUDA_UNSUPPORTED_SCALARS}|int" "SOUT=.*" "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
  # do int separately, because we will instantiate vector in it
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET_INT "S=int"                                         "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
  #
  ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_DOUBLE)
  IF(NOT KokkosClassic_ENABLE_CUDA_DOUBLE)
    APPEND_SET(Tpetra_ETI_EXCLUDE_SET "S=double SIN=double SOUT=double LO=.* GO=.* N=Kokkos::ThrustGPUNode")
  ENDIF()
  #
  ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_FLOAT)
  IF(NOT KokkosClassic_ENABLE_CUDA_FLOAT)
    APPEND_SET(Tpetra_ETI_EXCLUDE_SET "S=float SIN=float SOUT=float LO=.* GO=.* N=Kokkos::ThrustGPUNode")
  ENDIF()
ENDIF()

ASSERT_DEFINED(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
IF(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
  MESSAGE(STATUS "User/Downstream ETI set: ${Tpetra_ETI_LIBRARYSET}")
  MESSAGE(STATUS "Excluded instantiations: ${Tpetra_ETI_EXCLUDE_SET}:${Tpetra_ETI_EXCLUDE_SET_INT}")
  MESSAGE(STATUS "Full coverage explicit instantiation for scalars:         ${Tpetra_ETI_SCALARS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for global ordinals: ${Tpetra_ETI_GORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for local ordinals:  ${Tpetra_ETI_LORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for nodes:           ${Tpetra_ETI_NODES}")
  # generate ETI macros for Tpetra usage 
  # this uses the following lists: 
  #             nodes: Tpetra_ETI_NODES 
  #           scalars: Tpetra_ETI_SCALARS
  #   global ordinals: Tpetra_ETI_GORDS
  #    local ordinals: Tpetra_ETI_LORDS
  # assemble dual scalar (only) instantiations
  JOIN(Tpetra_ETI_LORDS   "|" FALSE ${Tpetra_ETI_LORDS}  )
  JOIN(Tpetra_ETI_GORDS   "|" FALSE ${Tpetra_ETI_GORDS}  )
  JOIN(Tpetra_ETI_NODES   "|" FALSE ${Tpetra_ETI_NODES}  )
  JOIN(Tpetra_ETI_SCALARS "|" FALSE ${Tpetra_ETI_SCALARS})
  # assemble single scalar instantiations
  TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts   "S=${Tpetra_ETI_SCALARS}" "N=${Tpetra_ETI_NODES}"
                                                 "LO=${Tpetra_ETI_LORDS}" "GO=${Tpetra_ETI_GORDS}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${SingleScalarInsts})
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "ETI set (before exclusions): ${Tpetra_ETI_LIBRARYSET}")
  ENDIF()
ELSEIF()
  # no ETI: these macros are used only for testing
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_LIBRARYSET "S=double" 
                                                   "LO=int" 
                                                   "GO=int" 
                                                   "N=${Tpetra_ETI_NODES}")
  IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
    TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_LIBRARYSET "S=std::complex<double>" 
                                                     "LO=int" 
                                                     "GO=int" 
                                                     "N=${Tpetra_ETI_NODES}")
  ENDIF()
ENDIF()

TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" 
    "${Tpetra_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_VECTOR(S,LO,GO,N)"            TPETRA_ETIMACRO_VECTOR)
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" 
    "${Tpetra_ETI_EXCLUDE_SET};${Tpetra_ETI_EXCLUDE_SET_INT}"  
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_SLGN(S,LO,GO,N)"            TPETRA_ETIMACRO_SLGN 
    "TPETRA_INSTANTIATE_LGN(LO,GO,N)"               TPETRA_ETIMACRO_LGN
    "TPETRA_INSTANTIATE_SLG(S,LO,GO)"               TPETRA_ETIMACRO_SLG  
    "TPETRA_INSTANTIATE_LG(LO,GO)"                  TPETRA_ETIMACRO_LG 
    "TPETRA_INSTANTIATE_N(N)"                       TPETRA_ETIMACRO_N
    "TPETRA_INSTANTIATE_TSLGN(CS,DS,LO,GO,N)"       TPETRA_ETIMACRO_TSLGN 
    "TPETRA_INSTANTIATE_TSLG(CS,DS,LO,GO)"          TPETRA_ETIMACRO_TSLG
    "TPETRA_INSTANTIATE_CONVERT(SOUT,SIN,LO,GO,N)"  TPETRA_ETIMACRO_CONVERT)
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" 
    "${Tpetra_ETI_EXCLUDE_SET};SIN=.* SOUT=.* CS=.* DS=.* S=.* LO=.* GO=.* N=Kokkos::ThrustGPUNode"  
    list_of_manglings   eti_typedefs
    "TPETRA_INSTANTIATE_SLGN_NOGPU(S,LO,GO,N)"            TPETRA_ETIMACRO_SLGN_NOGPU
    "TPETRA_INSTANTIATE_LGN_NOGPU(LO,GO,N)"               TPETRA_ETIMACRO_LGN_NOGPU
    "TPETRA_INSTANTIATE_SLG_NOGPU(S,LO,GO)"               TPETRA_ETIMACRO_SLG_NOGPU
    "TPETRA_INSTANTIATE_LG_NOGPU(LO,GO)"                  TPETRA_ETIMACRO_LG_NOGPU
    "TPETRA_INSTANTIATE_N_NOGPU(N)"                       TPETRA_ETIMACRO_N_NOGPU
    "TPETRA_INSTANTIATE_TSLGN_NOGPU(CS,DS,LO,GO,N)"       TPETRA_ETIMACRO_TSLGN_NOGPU
    "TPETRA_INSTANTIATE_TSLG_NOGPU(CS,DS,LO,GO)"          TPETRA_ETIMACRO_TSLG_NOGPU
    "TPETRA_INSTANTIATE_CONVERT_NOGPU(SOUT,SIN,LO,GO,N)"  TPETRA_ETIMACRO_CONVERT_NOGPU)
STRING(REPLACE "S="  "P=" ScalarToPacketSet "${Tpetra_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO=([^{ ]+)" "GO=\\1 P=\\1" GlobalToPacketSet1 "${Tpetra_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO={([^}]+)}" "GO={\\1} P={\\1}" GlobalToPacketSet2 "${Tpetra_ETI_LIBRARYSET}")
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}|P" "${ScalarToPacketSet};${GlobalToPacketSet1};${GlobalToPacketSet2}" 
    "${Tpetra_ETI_EXCLUDE_SET}"  
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_PLGN(P,LO,GO,N)"            TPETRA_ETIMACRO_PLGN)

# testing macros
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" 
    "${Tpetra_ETI_EXCLUDE_SET};S=int LO=.* GO=.* N=.*;S=long LO=.* GO=.* N=.*"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_TESTMV(S,LO,GO,N)"            TPETRA_ETIMACRO_TESTMV)
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" 
    "${Tpetra_ETI_EXCLUDE_SET};S=int LO=.* GO=.* N=.*;S=long LO=.* GO=.* N=.*;S=.* LO=.* GO=.* N=Kokkos::ThrustGPUNode"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_TESTMV_NOGPU(S,LO,GO,N)"      TPETRA_ETIMACRO_TESTMV_NOGPU)

TRIBITS_ETI_TYPE_EXPANSION(Tpetra_DII   "S=double" "N=${Tpetra_ETI_NODES}" "LO=int" "GO=int")
TRIBITS_ETI_GENERATE_MACROS(
    "${Tpetra_ETI_FIELDS}" "${Tpetra_DII}" 
    "${Tpetra_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(S,LO,GO,N)"      TPETRA_ETIMACRO_DII_NODE)

TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(TPETRA_ETI_TYPEDEFS "TPETRA_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

CONFIGURE_FILE(${Tpetra_SOURCE_DIR}/cmake/Tpetra_ETIHelperMacros.h.in ${Tpetra_BINARY_DIR}/src/Tpetra_ETIHelperMacros.h)
