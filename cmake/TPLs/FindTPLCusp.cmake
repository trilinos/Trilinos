INCLUDE(TPLDeclareLibraries)

IF (NOT TPL_ENABLE_Thrust)
  MESSAGE(FATAL_ERROR "\nCusp TPL requires that Thrust support is enabled. Please set \n  TPL_ENABLE_Thrust=ON\n\n")
ELSE()
  TPL_DECLARE_LIBRARIES( Cusp
    REQUIRED_HEADERS cusp/version.h
    )
ENDIF()
