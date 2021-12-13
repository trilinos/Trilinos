KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

# Workaround https://github.com/kokkos/kokkos/issues/4376 for ibm/xl
IF (NOT ${KOKKOS_COMPILER_IBM})
KOKKOSKERNELS_ADD_EXECUTABLE(
  gmres_ex_real_A
  SOURCES ex_real_A.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  gmres_test_real_A
  SOURCES test_real_A.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE(
  gmres_test_cmplx_A
  SOURCES test_cmplx_A.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  gmres_test_prec
  SOURCES test_prec.cpp
  )

ELSE ()
  MESSAGE (STATUS "SKIPPING gmres examples - Kokkos::complex<half_t> unsupported with ibm/xlC as host compiler")
ENDIF ()

