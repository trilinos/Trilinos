
# We need to inject the Trilinos/cmake directory to find
# TrilinosCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Trilinos_SOURCE_DIR}/cmake")

MACRO(TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11  CXX11_PACKAGE_NAME_IN)
  IF ("${${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}}" STREQUAL "")
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}=OFF"
      " because ${PROJECT_NAME}_ENABLE_CXX11='${${PROJECT_NAME}_ENABLE_CXX11}'!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN} OFF)
  ELSEIF (${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN})
    MESSAGE( FATAL_ERROR
      "ERROR: Setting"
      " ${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}='${${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}}'"
      " is not consistent with "
      " ${PROJECT_NAME}_ENABLE_CXX11='${${PROJECT_NAME}_ENABLE_CXX11}'!"
      " ${CXX11_PACKAGE_NAME_IN} requires C++11 support!  Either don't"
      " enable the package ${CXX11_PACKAGE_NAME_IN} or enable support for C++11!")
  ELSE()
    # This package is already disabled which is just fine.
  ENDIF()
ENDMACRO()


MACRO(TRIL_SET_BOOL_CACHE_VAR_FOR_CI  VAR_NAME  VAR_VALUE)
  IF ("${${VAR_NAME}}" STREQUAL "")
    MESSAGE("-- " "Setting ${VAR_NAME}='${VAR_VALUE}' by default for CI testing")
    SET(${VAR_NAME} ${VAR_VALUE} CACHE BOOL
      "Set in Trilinos/cmake/CallbackSetupExtraOptions.cmake")
  ENDIF()
ENDMACRO()


MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  #MESSAGE("TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS got called!")

  SET(TPL_ENABLE_MPI OFF CACHE BOOL "Enable MPI support.")

  #
  # Set options for global enable/disable of float and complex
  #

  SET(Trilinos_ENABLE_FLOAT  OFF  CACHE  BOOL
    "Enable the float scalar type in all Trilinos packages by default.")

  SET(Trilinos_ENABLE_COMPLEX  OFF  CACHE  BOOL
    "Enable std::complex<T> scalar types in all Trilinos packages by default.")

  IF (Trilinos_ENABLE_COMPLEX  AND  Trilinos_ENABLE_FLOAT)
    SET(Trilinos_ENABLE_COMPLEX_FLOAT_DEFAULT  ON)
  ELSE()
    SET(Trilinos_ENABLE_COMPLEX_FLOAT_DEFAULT  OFF)
  ENDIF()
  SET(Trilinos_ENABLE_COMPLEX_FLOAT  ${Trilinos_ENABLE_COMPLEX_FLOAT_DEFAULT}
    CACHE  BOOL
    "Enable std::complex<float> scalar types in all Trilinos packages by default.")

  SET(Trilinos_ENABLE_COMPLEX_DOUBLE  ${Trilinos_ENABLE_COMPLEX}
    CACHE  BOOL
    "Enable std::complex<double> scalar types in all Trilinos packages by default.")

  #
  # Trilinos Data Dir?  Is this still being used anywhere?
  #

  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )

  #
  # Put in disables based on various criteria
  #
    
  IF (NOT ${PROJECT_NAME}_ENABLE_CXX11)
    TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11("Kokkos")
    TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11("Tpetra")
  ENDIF()
    
  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_ForTrilinos=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_ForTrilinos OFF)
  ENDIF()

  IF ("${${PROJECT_NAME}_ENABLE_PyTrilinos}" STREQUAL "" AND NOT BUILD_SHARED_LIBS)
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_PyTrilinos=OFF"
      " because BUILD_SHARED_LIBS=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF)
  ENDIF()

  IF (NOT EXISTS "${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota")
    MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_TriKota=OFF"
      " because '${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota' does not exist!")
    SET(${PROJECT_NAME}_ENABLE_TriKota OFF)
  ENDIF()
    
  # Used by some Trilinos packages?
  SET(TRILINOS_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

  #
  # Adjust options for a standard CI build
  #

  ADVANCED_SET(Trilinos_ENABLE_CI_TEST_MODE  OFF
    CACHE  BOOL
    "If set to 'ON', then Trilinos packages will be put into CI mode which adjusts a bunch fo things"
    )

  IF (Trilinos_ENABLE_CI_TEST_MODE)

    MESSAGE(
      "-- "
      "NOTE: Trilinos_ENABLE_CI_TEST_MODE='${Trilinos_ENABLE_CI_TEST_MODE}'"
      ", setting cache vars for Trilinos CI testing mode ..."
      )

    # Turn off float and complex by default
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Teuchos_ENABLE_FLOAT OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Teuchos_ENABLE_COMPLEX OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Sacado_ENABLE_COMPLEX OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Thyra_ENABLE_COMPLEX OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Tpetra_INST_COMPLEX_DOUBLE OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Tpetra_INST_COMPLEX_FLOAT OFF)
    #TRIL_SET_BOOL_CACHE_VAR_FOR_CI(Anasazi_ENABLE_COMPLEX OFF)
    # ToDo: Remove the above once Trlinos_ENABLE_FLOAT and
    # Trilinos_ENABLE_COMPLEX are supported and are off by default (see
    # Trilinos GitHub #362)

    # Disable a bunch of TPLs that are not to be enabled in CI testing (since
    # the SEMS env does not have them).
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_GLM OFF)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_Matio OFF)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_X11 OFF)

    # Default enable the TPLs that SEMS provides
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_Pthread ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_BLAS ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_LAPACK ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_Boost ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_BoostLib ON)
    IF (TPL_ENABLE_MPI)
      TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_ParMETIS ON)
    ENDIF()
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_Zlib ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_HDF5 ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_Netcdf ON)
    TRIL_SET_BOOL_CACHE_VAR_FOR_CI(TPL_ENABLE_SuperLU ON)

    # Disable long-failing Pir test until it can be fixed (#826)
    SET(Piro_EpetraSolver_MPI_4_DISABLE ON)

  ENDIF()

  # NOTE: Above, the cache var Trilinos_ENABLE_CI_TEST_MODE and the above code
  # helps to remove dupliction between the checkin-test.py script and the
  # post-push CI server.  This could be moved into a *.cmake fragment file
  # however.

ENDMACRO()
