
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

  OPTION(Trilinos_ENABLE_THREAD_SAFE
    "Enable thread safe code including RCP classes." OFF )

  ASSERT_DEFINED(${PROJECT_NAME}_ENABLE_CXX11)
  IF (Trilinos_ENABLE_THREAD_SAFE AND NOT ${PROJECT_NAME}_ENABLE_CXX11)
    MESSAGE(FATAL_ERROR
      "You set Trilinos_ENABLE_THREAD_SAFE=ON, but ${PROJECT_NAME}' support"
      " for CXX11 is not enabled (${PROJECT_NAME}_ENABLE_CXX11=OFF)."
      "  This is not allowed.  Please enable ${PROJECT_NAME}_ENABLE_CXX11 in"
      " ${PROJECT_NAME} before attempting to enable Trilinos_ENABLE_THREAD_SAFE"
      " or leave Trilinos_ENABLE_THREAD_SAFE off.")
  ENDIF ()

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
    
  IF (
      NOT ${PROJECT_NAME}_ENABLE_Fortran
      AND
      (
        "${${PROJECT_NAME}_ENABLE_ForTrilinos}" STREQUAL ""
        OR
        ${PROJECT_NAME}_ENABLE_ForTrilinos
      )
    )
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_ForTrilinos=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_ForTrilinos OFF)
  ENDIF()

  IF (
      NOT BUILD_SHARED_LIBS
      AND
      (
        "${${PROJECT_NAME}_ENABLE_PyTrilinos}" STREQUAL ""
        OR
        ${PROJECT_NAME}_ENABLE_PyTrilinos
      )
    )
    MESSAGE(
      "\n***"
      "\n*** NOTE: Setting ${PROJECT_NAME}_ENABLE_PyTrilinos=OFF"
      " because BUILD_SHARED_LIBS=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF)
  ENDIF()

  IF (
      NOT EXISTS "${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota"
      AND
      (
        "${${PROJECT_NAME}_ENABLE_TriKota}" STREQUAL ""
        OR
        ${PROJECT_NAME}_ENABLE_TriKota
      )
    )
    MESSAGE("-- " "Setting ${PROJECT_NAME}_ENABLE_TriKota=OFF"
      " because '${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota' does not exist!")
    SET(${PROJECT_NAME}_ENABLE_TriKota OFF)
  ENDIF()
    
  # Used by some Trilinos packages?
  SET(TRILINOS_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

ENDMACRO()
