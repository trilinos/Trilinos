# Tpetra ETI type fields
SET(Tpetra_ETI_FIELDS "CS;DS;LO;GO;N" CACHE STRING "")

IF(Tpetra_ENABLE_EXPLICIT_INSTANTIATION)
  IF(Tpetra_ENABLE_Thrust)
    TpetraExpandTypesetProduct(Tpetra_ETI_EXCLUDE_SET "qd_real;dd_real;int" ".*" ".*" "Kokkos::ThrustGPUNode")
    IF(NOT KokkosClassic_ENABLE_CUDA_DOUBLE)
      TpetraExpandTypesetProduct(Tpetra_ETI_EXCLUDE_SET "double" ".*" ".*" "Kokkos::ThrustGPUNode")
    ENDIF()
    IF(NOT KokkosClassic_ENABLE_CUDA_FLOAT)
      TpetraExpandTypesetProduct(Tpetra_ETI_EXCLUDE_SET "float"  ".*" ".*" "Kokkos::ThrustGPUNode")
    ENDIF()
  ENDIF()

  MESSAGE(STATUS "User-specified ETI set: ${Tpetra_ETI_LIBRARYSET}")
  MESSAGE(STATUS "Excluded instantiations: ${Tpetra_ETI_EXCLUDE_SET}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following nodes: ${Tpetra_ETI_NODES}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following scalars: ${Tpetra_ETI_SCALARS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following global ordinals: ${Tpetra_ETI_GORDS}")
  MESSAGE(STATUS "Full coverage explicit instantiation for following local ordinals: ${Tpetra_ETI_LORDS}")
  # generate ETI macros for Tpetra usage 
  # this uses the following lists: 
  #             nodes: Tpetra_ETI_NODES 
  #           scalars: Tpetra_ETI_SCALARS
  #   global ordinals: Tpetra_ETI_GORDS
  #    local ordinals: Tpetra_ETI_LORDS
  TpetraExpandTypesetProduct(Tpetra_ETI_LIBRARYSET "${Tpetra_ETI_SCALARS}" "${Tpetra_ETI_LORDS}" "${Tpetra_ETI_GORDS}" "${Tpetra_ETI_NODES}")
  MESSAGE(STATUS "ETI set before exclusions: ${Tpetra_ETI_LIBRARYSET}")
ELSEIF()
  # these macros are used only for testing
  TpetraExpandTypesetProduct(Tpetra_ETI_LIBRARYSET "int;double" "int" "int;long" "${Tpetra_ETI_NODES}")
  IF (SS_FOR_DEV_PS_FOR_RELEASE AND HAVE_COMPLEX_BLAS)
    TpetraExpandTypesetProduct(Tpetra_ETI_LIBRARYSET "std::complex<double>" "int" "int;long" "${Tpetra_ETI_NODES}")
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
