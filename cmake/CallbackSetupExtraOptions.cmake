
# We need to inject the Trilinos/cmake directory to find
# TrilinosCreateClientTemplateHeaders.cmake
SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} "${Trilinos_SOURCE_DIR}/cmake")

MACRO(TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS)

  #MESSAGE("TRIBITS_REPOSITORY_SETUP_EXTRA_OPTIONS got called!")

  SET(TPL_ENABLE_MPI OFF CACHE BOOL "Enable MPI support.")

  #
  # Set options for global enable/disable of float and complex
  #

  SET(Trilinos_ENABLE_FLOAT  OFF  CACHE  BOOL
    "Enable the float scalar type in all Trilinos packages by default.")

  SET(Trilinos_ENABLE_LONG_DOUBLE  OFF  CACHE  BOOL
    "Enable the long double scalar type in all Trilinos packages by default.")

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

  #
  # Trilinos Data Dir?  Is this still being used anywhere?
  #

  ADVANCED_SET(Trilinos_DATA_DIR  NOTFOUND
    CACHE PATH
    "Path TrilinosData directory to find more tests and other stuff" )

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


#
# Set up for build stats
#

include("${Trilinos_SOURCE_DIR}/commonTools/build_stats/BuildStatsWrappers.cmake")
generate_build_stats_wrappers()
remove_build_stats_file_on_configure()
remove_build_stats_timing_files_on_fresh_configure()

# Set up C++ language standard selection
include(SetTrilinosCxxStandard)
