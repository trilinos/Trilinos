# Tpetra ETI type fields
SET(Tpetra_ETI_FIELDS "CS|DS|LO|GO|N" CACHE STRING "")

ASSERT_DEFINED(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
IF(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
  ASSERT_DEFINED(Tpetra_ENABLE_Thrust)
  IF(Tpetra_ENABLE_Thrust)
    TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET "CS=qd_real|dd_real|int" "DS=qd_real|dd_real|int"
                                                      "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
    #
    ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_DOUBLE)
    IF(NOT KokkosClassic_ENABLE_CUDA_DOUBLE)
      TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET "CS=double" "DS=double"
                                                        "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
    ENDIF()
    #
    ASSERT_DEFINED(KokkosClassic_ENABLE_CUDA_FLOAT)
    IF(NOT KokkosClassic_ENABLE_CUDA_FLOAT)
      TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_EXCLUDE_SET "CS=float" "DS=float"
                                                        "LO=.*" "GO=.*" "N=Kokkos::ThrustGPUNode")
    ENDIF()
  ENDIF()

  MESSAGE(STATUS "User/Downstream ETI set: ${Tpetra_ETI_LIBRARYSET}")
  MESSAGE(STATUS "Excluded instantiations: ${Tpetra_ETI_EXCLUDE_SET}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following scalars:         ${Tpetra_ETI_SCALARS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following global ordinals: ${Tpetra_ETI_GORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following local ordinals:  ${Tpetra_ETI_LORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following nodes:           ${Tpetra_ETI_NODES}")
  # generate ETI macros for Tpetra usage 
  # this uses the following lists: 
  #             nodes: Tpetra_ETI_NODES 
  #           scalars: Tpetra_ETI_SCALARS
  #   global ordinals: Tpetra_ETI_GORDS
  #    local ordinals: Tpetra_ETI_LORDS
  JOIN(Tpetra_ETI_SCALARS "|" FALSE ${Tpetra_ETI_SCALARS})
  JOIN(Tpetra_ETI_LORDS   "|" FALSE ${Tpetra_ETI_LORDS}  )
  JOIN(Tpetra_ETI_GORDS   "|" FALSE ${Tpetra_ETI_GORDS}  )
  JOIN(Tpetra_ETI_NODES   "|" FALSE ${Tpetra_ETI_NODES}  )
  TRIBITS_ETI_TYPE_EXPANSION(FullCoverageInsts "CS=${Tpetra_ETI_SCALARS}" "DS=${Tpetra_ETI_SCALARS}" 
                                          "LO=${Tpetra_ETI_LORDS}" "GO=${Tpetra_ETI_GORDS}" 
                                          "N=${Tpetra_ETI_NODES}")
  TRIBITS_ADD_ETI_INSTANTIATIONS(Tpetra ${FullCoverageInsts})
  MESSAGE(STATUS "Instantiation set before exclusions: ${Tpetra_ETI_LIBRARYSET}")
ELSEIF()
  # these macros are used only for testing
  TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_LIBRARYSET "CS=int|double" "DS=int|double" 
                                                   "LO=int" "GO=int|long" "N=${Tpetra_ETI_NODES}")
  IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
    TRIBITS_ETI_TYPE_EXPANSION(Tpetra_ETI_LIBRARYSET "CS=std::complex<double>" "DS=std::complex<double>" 
                                                     "LO=int" "GO=int|long" "N=${Tpetra_ETI_NODES}")
  ENDIF()
ENDIF()

TRIBITS_GENERATE_ETI_MACROS("${Tpetra_ETI_FIELDS}" "${Tpetra_ETI_LIBRARYSET}" "${Tpetra_ETI_EXCLUDE_SET}"  
                            "TPETRA_INSTANTIATE_TSLGN(CS,DS,LO,GO,N)" TPETRA_ETIMACRO_TSLGN 
                            "TPETRA_INSTANTIATE_SLGN(DS,LO,GO,N)"     TPETRA_ETIMACRO_SLGN 
                            "TPETRA_INSTANTIATE_LGN(LO,GO,N)"         TPETRA_ETIMACRO_LGN
                            "TPETRA_INSTANTIATE_TSLG(CS,DS,LO,GO)"    TPETRA_ETIMACRO_TSLG  
                            "TPETRA_INSTANTIATE_SLG(DS,LO,GO)"        TPETRA_ETIMACRO_SLG  
                            "TPETRA_INSTANTIATE_LG(LO,GO)"            TPETRA_ETIMACRO_LG 
                            "TPETRA_INSTANTIATE_N(N)"                 TPETRA_ETIMACRO_N)

CONFIGURE_FILE(${Tpetra_SOURCE_DIR}/cmake/Tpetra_ETIHelperMacros.h.in ${Tpetra_BINARY_DIR}/src/Tpetra_ETIHelperMacros.h)
