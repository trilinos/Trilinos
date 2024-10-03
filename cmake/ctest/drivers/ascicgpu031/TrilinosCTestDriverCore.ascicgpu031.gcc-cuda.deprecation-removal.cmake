
INCLUDE("${CTEST_SCRIPT_DIRECTORY}/../../TrilinosCTestDriverCore.cmake")

#
# Platform/compiler specific options for geminga using gcc
#

MACRO(TRILINOS_SYSTEM_SPECIFIC_CTEST_DRIVER)

  # Base of Trilinos/cmake/ctest then BUILD_DIR_NAME
  IF(COMM_TYPE STREQUAL MPI)
    string(TOUPPER $ENV{SIERRA_MPI} UC_MPI_NAME)
    SET(BUILD_DIR_NAME ${UC_MPI_NAME}_${BUILD_TYPE}_${BUILD_NAME_DETAILS})
  ELSE()
    SET(BUILD_DIR_NAME ${COMM_TYPE}-${BUILD_TYPE}_${BUILD_NAME_DETAILS})
  ENDIF()

  SET(Trilinos_REPOSITORY_LOCATION_NIGHTLY_DEFAULT "git@gitlab-ex.sandia.gov:trilinos-project/Trilinos.git")
  SET(Trilinos_BRANCH "tpetra-deprecation-removal" )

  SET(CTEST_DASHBOARD_ROOT  "${TRILINOS_CMAKE_DIR}/../../${BUILD_DIR_NAME}" )
  SET(CTEST_NOTES_FILES     "${CTEST_SCRIPT_DIRECTORY}/${CTEST_SCRIPT_NAME}" )
  SET(CTEST_BUILD_FLAGS     "-j45 -i" )

  SET_DEFAULT(CTEST_PARALLEL_LEVEL                  "45" )
  SET_DEFAULT(Trilinos_ENABLE_SECONDARY_TESTED_CODE ON)
  SET(Trilinos_CTEST_DO_ALL_AT_ONCE FALSE)
  SET_DEFAULT(Trilinos_EXCLUDE_PACKAGES             ${EXTRA_EXCLUDE_PACKAGES} TriKota Optika Pamgen)

  # Select package disables
  set (Trilinos_ENABLE_Gtest OFF CACHE BOOL "Gtest just does not build" FORCE)
  set (Trilinos_ENABLE_ShyLU_NodeTacho OFF CACHE BOOL "Can't test Tacho with CUDA without RDC" FORCE)
  set (Trilinos_ENABLE_Shards OFF CACHE BOOL "Shards does not build" FORCE)
  set (Trilinos_ENABLE_Epetra OFF CACHE BOOL "We do not want Epetra" FORCE)

  # Select test disables
  set (Kokkos_CoreUnitTest_CudaTimingBased_MPI_1_DISABLE ON CACHE BOOL "Not to be run in nightly testing" FORCE)

  SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
      "-DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE}"

      # Adding the following as a possible fix for github issue #2115.
      #KDD This flag appears to be unnecessary in April 2021, and it
      #KDD breaks building of Zoltan tests
      #KDD "-DCMAKE_CXX_USE_RESPONSE_FILE_FOR_OBJECTS:BOOL=ON"

      ### ALWAYS AND EVERYWHERE ###
      "-DTrilinos_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON"
      "-DBUILD_SHARED_LIBS:BOOL=ON"
      "-DTrilinos_ENABLE_TESTS:BOOL=ON"
      "-DTrilinos_ENABLE_EXAMPLES:BOOL=ON"
      "-DTrilinos_ENABLE_DEPENDENCY_UNIT_TESTS:BOOL=OFF"
      "-DTeuchos_GLOBALLY_REDUCE_UNITTEST_RESULTS:BOOL=ON"
      "-DTrilinos_ENABLE_COMPLEX=ON"
      "-DTeuchos_ENABLE_COMPLEX=ON"
      "-DTpetra_INST_COMPLEX_DOUBLE=ON"

      ### COMPILERS AND FLAGS ###
      "-DCMAKE_CXX_FLAGS:STRING='-Wall -Wno-unknown-pragmas -Wno-unused-but-set-variable -Wno-inline -Wshadow'"
      "-DTrilinos_ENABLE_Fortran:BOOL=OFF"

      ### TPLS ###
      "-DTPL_ENABLE_CUDA:BOOL=ON"
      "-DCMAKE_POLICY_DEFAULT_CMP0074=NEW"
      "-DTPL_ENABLE_CUSPARSE:BOOL=ON"
      "-DTPL_ENABLE_HWLOC:BOOL=OFF"


      # Host Blas is required (https://github.com/kokkos/kokkos-kernels/issues/347) for Kokkos-Kernels to build correctly
      "-DTPL_ENABLE_BLAS:BOOL=ON"
      "-DTPL_ENABLE_LAPACK:BOOL=ON"
      "-DTPL_BLAS_LIBRARIES=/usr/lib64/libblas.so"
      "-DTPL_LAPACK_LIBRARIES=/usr/lib64/liblapack.so"

      ### PACKAGE CONFIGURATION ###
      "-DKokkos_ENABLE_CUDA:BOOL=ON"
      "-DKokkos_ENABLE_CUDA_LAMBDA:BOOL=ON"
      "-DKokkos_ARCH_SKX:BOOL=ON"
      "-DKokkos_ARCH_VOLTA70:BOOL=ON"
      "-DKokkos_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE=OFF"

      "-DTrilinos_ENABLE_Epetra:BOOL=OFF"
      "-DTrilinos_ENABLE_Gtest:BOOL=OFF"
      "-DTrilinos_ENABLE_Pamgen:BOOL=OFF"
      "-DTrilinos_ENABLE_Shards:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_Node:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_NodeTacho:BOOL=OFF"            
      "-DTrilinos_ENABLE_ShyLU:BOOL=OFF"
      "-DTrilinos_ENABLE_ShyLU_DD:BOOL=OFF"
      "-DAmesos2_ENABLE_ShyLU_NodeTacho:BOOL=OFF"
      "-DAmesos2_ENABLE_ShyLU_NodeBasker:BOOL=OFF"

      ### MISC ###
      "-DCMAKE_VERBOSE_MAKEFILE:BOOL=ON"
  )

  SET_DEFAULT(COMPILER_VERSION "$ENV{SIERRA_PLATFORM}")

  # Ensure that MPI is on for all parallel builds that might be run.
  IF(COMM_TYPE STREQUAL MPI)

    SET(EXTRA_SYSTEM_CONFIGURE_OPTIONS
        ${EXTRA_SYSTEM_CONFIGURE_OPTIONS}
        "-DTPL_ENABLE_MPI:BOOL=ON"
        "-DMPI_BASE_DIR:PATH=$ENV{MPIHOME}"
       )

  ENDIF()

  TRILINOS_CTEST_DRIVER()

ENDMACRO()
