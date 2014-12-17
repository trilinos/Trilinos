# The file that is produced by this module
SET(${PACKAGE_NAME}_ETI_FILE TpetraCore_ETIHelperMacros.h)
SET(${PACKAGE_NAME}_ETI_FILE_PATH ${Tpetra_BINARY_DIR}/core/src/${${PACKAGE_NAME}_ETI_FILE})

#
# A) See if a static pre-created file is provide and us it if it is
#

ADVANCED_SET(Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE ""
  CACHE PATH
  "If set, gives the path to a static version of the file ${${PACKAGE_NAME}_ETI_FILE}.  If not set (default '') then the file is generated automatically (and at great cost)"
  )

IF(Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE)
  MESSAGE("-- NOTE: Skipping generation and using provided static file"
     " '${Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE}'")
  CONFIGURE_FILE(
    ${Tpetra_USE_STATIC_ETI_MACROS_HEADER_FILE}
    ${${PACKAGE_NAME}_ETI_FILE_PATH}
    COPYONY
    )
  RETURN()
ENDIF()

#
# B) We must generate the file anew :-(
#

# Tpetra ETI type fields
SET(${PACKAGE_NAME}_ETI_FIELDS "SIN|SOUT|S|LO|GO|N|CS|DS")

# Exclude all of the types that CUDA/Thrust doesn't support
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Thrust)
IF(${PACKAGE_NAME}_ENABLE_Thrust)
  # no dd_real/qd_real support for CUDA, nor int/complex even via Cusp :( 
  SET(CUDA_UNSUPPORTED_SCALARS "dd_real|qd_real|std::complex<double>|std::complex<float>")
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET     "S=${CUDA_UNSUPPORTED_SCALARS}"                 "LO=.*" "GO=.*" "N=KokkosClassic::ThrustGPUNode")
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET     "SIN=.*" "SOUT=${CUDA_UNSUPPORTED_SCALARS}|int" "LO=.*" "GO=.*" "N=KokkosClassic::ThrustGPUNode")
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET     "SIN=${CUDA_UNSUPPORTED_SCALARS}|int" "SOUT=.*" "LO=.*" "GO=.*" "N=KokkosClassic::ThrustGPUNode")
  # do int separately, because we will instantiate vector in it
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_EXCLUDE_SET_INT "S=int|long|unsigned int"                       "LO=.*" "GO=.*" "N=KokkosClassic::ThrustGPUNode")
  #
  ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_DOUBLE)
  IF(NOT KokkosClassic_ENABLE_CUDA_DOUBLE)
    APPEND_SET(${PACKAGE_NAME}_ETI_EXCLUDE_SET "S=double SIN=double SOUT=double LO=.* GO=.* N=KokkosClassic::ThrustGPUNode")
  ENDIF()
  #
  ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_FLOAT)
  IF(NOT KokkosClassic_ENABLE_CUDA_FLOAT)
    APPEND_SET(${PACKAGE_NAME}_ETI_EXCLUDE_SET "S=float SIN=float SOUT=float LO=.* GO=.* N=KokkosClassic::ThrustGPUNode")
  ENDIF()
ENDIF()

ASSERT_DEFINED(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
IF(TpetraCore_ENABLE_EXPLICIT_INSTANTIATION)
  MESSAGE(STATUS "User/Downstream ETI set: ${${PACKAGE_NAME}_ETI_LIBRARYSET}")
  MESSAGE(STATUS "Excluded instantiations: ${${PACKAGE_NAME}_ETI_EXCLUDE_SET}:${${PACKAGE_NAME}_ETI_EXCLUDE_SET_INT}")
  MESSAGE(STATUS "Full coverage explicit instantiation for scalars:         ${${PACKAGE_NAME}_ETI_SCALARS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for global ordinals: ${${PACKAGE_NAME}_ETI_GORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for local ordinals:  ${${PACKAGE_NAME}_ETI_LORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for nodes:           ${${PACKAGE_NAME}_ETI_NODES}")
  # generate ETI macros for Tpetra usage 
  # this uses the following lists: 
  #             nodes: ${PACKAGE_NAME}_ETI_NODES 
  #           scalars: ${PACKAGE_NAME}_ETI_SCALARS
  #   global ordinals: ${PACKAGE_NAME}_ETI_GORDS
  #    local ordinals: ${PACKAGE_NAME}_ETI_LORDS
  # assemble dual scalar (only) instantiations
  JOIN(${PACKAGE_NAME}_ETI_LORDS   "|" FALSE ${${PACKAGE_NAME}_ETI_LORDS}  )
  JOIN(${PACKAGE_NAME}_ETI_GORDS   "|" FALSE ${${PACKAGE_NAME}_ETI_GORDS}  )
  JOIN(${PACKAGE_NAME}_ETI_NODES   "|" FALSE ${${PACKAGE_NAME}_ETI_NODES}  )
  JOIN(${PACKAGE_NAME}_ETI_SCALARS "|" FALSE ${${PACKAGE_NAME}_ETI_SCALARS})
  # assemble single scalar instantiations
  TRIBITS_ETI_TYPE_EXPANSION(SingleScalarInsts   "S=${${PACKAGE_NAME}_ETI_SCALARS}" "N=${${PACKAGE_NAME}_ETI_NODES}"
                                                 "LO=${${PACKAGE_NAME}_ETI_LORDS}" "GO=${${PACKAGE_NAME}_ETI_GORDS}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(TpetraCore ${SingleScalarInsts})
  IF(${PROJECT_NAME}_VERBOSE_CONFIGURE)
    MESSAGE(STATUS "ETI set (before exclusions): ${${PACKAGE_NAME}_ETI_LIBRARYSET}")
  ENDIF()
ELSE()
  # no ETI: these macros are used only for testing
  JOIN(${PACKAGE_NAME}_ETI_NODES   "|" FALSE ${${PACKAGE_NAME}_ETI_NODES}  )
  TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_LIBRARYSET "S=double" 
                                                   "LO=int" 
                                                   "GO=int" 
                                                   "N=${${PACKAGE_NAME}_ETI_NODES}")
  IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
    TRIBITS_ETI_TYPE_EXPANSION(${PACKAGE_NAME}_ETI_LIBRARYSET "S=std::complex<double>" 
                                                     "LO=int" 
                                                     "GO=int" 
                                                     "N=${${PACKAGE_NAME}_ETI_NODES}")
  ENDIF()
ENDIF()

TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_VECTOR(S,LO,GO,N)"            TPETRA_ETIMACRO_VECTOR)
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};${${PACKAGE_NAME}_ETI_EXCLUDE_SET_INT}"  
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_SLGN(S,LO,GO,N)"            TPETRA_ETIMACRO_SLGN 
    "TPETRA_INSTANTIATE_LGN(LO,GO,N)"               TPETRA_ETIMACRO_LGN
    "TPETRA_INSTANTIATE_SLG(S,LO,GO)"               TPETRA_ETIMACRO_SLG  
    "TPETRA_INSTANTIATE_LG(LO,GO)"                  TPETRA_ETIMACRO_LG 
    "TPETRA_INSTANTIATE_SL(S,LO)"                   TPETRA_ETIMACRO_SL
    "TPETRA_INSTANTIATE_N(N)"                       TPETRA_ETIMACRO_N
    "TPETRA_INSTANTIATE_TSLGN(CS,DS,LO,GO,N)"       TPETRA_ETIMACRO_TSLGN 
    "TPETRA_INSTANTIATE_TSLG(CS,DS,LO,GO)"          TPETRA_ETIMACRO_TSLG
    "TPETRA_INSTANTIATE_CONVERT(SOUT,SIN,LO,GO,N)"  TPETRA_ETIMACRO_CONVERT)
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};SIN=.* SOUT=.* CS=.* DS=.* S=.* LO=.* GO=.* N=KokkosClassic::ThrustGPUNode"  
    list_of_manglings   eti_typedefs
    "TPETRA_INSTANTIATE_SLGN_NOGPU(S,LO,GO,N)"            TPETRA_ETIMACRO_SLGN_NOGPU
    "TPETRA_INSTANTIATE_LGN_NOGPU(LO,GO,N)"               TPETRA_ETIMACRO_LGN_NOGPU
    "TPETRA_INSTANTIATE_SLG_NOGPU(S,LO,GO)"               TPETRA_ETIMACRO_SLG_NOGPU
    "TPETRA_INSTANTIATE_LG_NOGPU(LO,GO)"                  TPETRA_ETIMACRO_LG_NOGPU
    "TPETRA_INSTANTIATE_N_NOGPU(N)"                       TPETRA_ETIMACRO_N_NOGPU
    "TPETRA_INSTANTIATE_TSLGN_NOGPU(CS,DS,LO,GO,N)"       TPETRA_ETIMACRO_TSLGN_NOGPU
    "TPETRA_INSTANTIATE_TSLG_NOGPU(CS,DS,LO,GO)"          TPETRA_ETIMACRO_TSLG_NOGPU
    "TPETRA_INSTANTIATE_CONVERT_NOGPU(SOUT,SIN,LO,GO,N)"  TPETRA_ETIMACRO_CONVERT_NOGPU)
STRING(REPLACE "S="  "P=" ScalarToPacketSet "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO=([^{ ]+)" "GO=\\1 P=\\1" GlobalToPacketSet1 "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
STRING(REGEX REPLACE "GO={([^}]+)}" "GO={\\1} P={\\1}" GlobalToPacketSet2 "${${PACKAGE_NAME}_ETI_LIBRARYSET}")
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}|P" "${ScalarToPacketSet};${GlobalToPacketSet1};${GlobalToPacketSet2}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"  
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_PLGN(P,LO,GO,N)"            TPETRA_ETIMACRO_PLGN)

# testing macros
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};S=int LO=.* GO=.* N=.*;S=long LO=.* GO=.* N=.*"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_TESTMV(S,LO,GO,N)"            TPETRA_ETIMACRO_TESTMV)
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${${PACKAGE_NAME}_ETI_LIBRARYSET}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET};S=int LO=.* GO=.* N=.*;S=long LO=.* GO=.* N=.*;S=.* LO=.* GO=.* N=KokkosClassic::ThrustGPUNode"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_TESTMV_NOGPU(S,LO,GO,N)"      TPETRA_ETIMACRO_TESTMV_NOGPU)

TRIBITS_ETI_TYPE_EXPANSION(Tpetra_DII   "S=double" "N=${${PACKAGE_NAME}_ETI_NODES}" "LO=int" "GO=int")
TRIBITS_ETI_GENERATE_MACROS(
    "${${PACKAGE_NAME}_ETI_FIELDS}" "${Tpetra_DII}" 
    "${${PACKAGE_NAME}_ETI_EXCLUDE_SET}"
    list_of_manglings eti_typedefs
    "TPETRA_INSTANTIATE_DOUBLE_INT_INT_N(S,LO,GO,N)"      TPETRA_ETIMACRO_DII_NODE)

TRIBITS_ETI_GENERATE_TYPEDEF_MACRO(TPETRA_ETI_TYPEDEFS "TPETRA_ETI_MANGLING_TYPEDEFS" "${eti_typedefs}")

CONFIGURE_FILE(
  ${${PACKAGE_NAME}_SOURCE_DIR}/cmake/TpetraCore_ETIHelperMacros.h.in
  ${${PACKAGE_NAME}_ETI_FILE_PATH}
  )
