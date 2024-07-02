
IF (NOT TPL_ENABLE_Thrust)
  MESSAGE(FATAL_ERROR "\nCusp TPL requires that Thrust support is enabled. Please set \n  TPL_ENABLE_Thrust=ON\n\n")
ELSE()
  TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( Cusp
    REQUIRED_HEADERS cusp/complex.h cusp/csr_matrix.h
    )
ENDIF()
