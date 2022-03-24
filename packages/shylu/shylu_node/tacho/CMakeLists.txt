TRIBITS_SUBPACKAGE(Tacho)

IF (Kokkos_ENABLE_CUDA)
  # Tacho requires CUDA >= 8.0.
  IF (DEFINED CUDA_VERSION AND (CUDA_VERSION VERSION_LESS "8.0"))
    MESSAGE(FATAL_ERROR "Tacho requires CUDA 8 if CUDA is enabled")
  ENDIF()
ENDIF()

IF (Kokkos_ENABLE_PTHREAD)
  IF (NOT Kokkos_ENABLE_OPENMP)
    MESSAGE(FATAL_ERROR "Tacho can not be build with Pthreads as the Kokkos Host Backend.")
  ENDIF()
ENDIF()

ADD_SUBDIRECTORY(src)
IF (Trilinos_ENABLE_Gtest)
  TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)
  #TRIBITS_ADD_TEST_DIRECTORIES(unit-test)
ELSE()
  MESSAGE(STATUS "Tacho disables examples and tests as Trilinos disables Gtest")
ENDIF()

TRIBITS_SUBPACKAGE_POSTPROCESS()
