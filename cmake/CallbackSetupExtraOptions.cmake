
# We need to inject the Trilinos/cmake directory to find
# TrilinosCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Trilinos_SOURCE_DIR}/cmake")

# Install TriBITS so that other projects can use it.
INSTALL(
  DIRECTORY "${Trilinos_SOURCE_DIR}/cmake/tribits"
  DESTINATION "${${CMAKE_PROJECT_NAME}_INSTALL_LIB_DIR}/cmake"
  )

MACRO(TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11  CXX11_PACKAGE_NAME_IN)
  IF ("${${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}}" STREQUAL "")
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_${CXX11_PACKAGE_NAME_IN}=OFF"
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

  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )
    
  IF (NOT ${PROJECT_NAME}_ENABLE_CXX11)
    TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11("Kokkos")
    TRILINOS_DISABLE_PACKAGE_REQUIRING_CXX11("Tpetra")
  ENDIF()
    
  IF (NOT ${PROJECT_NAME}_ENABLE_Fortran)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_ForTrilinos=OFF"
      " because ${PROJECT_NAME}_ENABLE_Fortran=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_ForTrilinos OFF)
  ENDIF()

  IF ("${${PROJECT_NAME}_ENABLE_PyTrilinos}" STREQUAL "" AND NOT BUILD_SHARED_LIBS)
    MESSAGE(
      "\n***"
      "\n*** Warning: Setting ${PROJECT_NAME}_ENABLE_PyTrilinos=OFF"
      " because BUILD_SHARED_LIBS=OFF!"
      "\n***\n"
      )
    SET(${PROJECT_NAME}_ENABLE_PyTrilinos OFF)
  ENDIF()

  IF (NOT EXISTS "${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota")
    MESSAGE("-- " "  Setting ${PROJECT_NAME}_ENABLE_TriKota=OFF"
      " because '${Trilinos_SOURCE_DIR}/packages/TriKota/Dakota' does not exit!")
    SET(${PROJECT_NAME}_ENABLE_TriKota OFF)
  ENDIF()
    
  # Used by some Trilinos packages?
  SET(TRILINOS_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})

ENDMACRO()
