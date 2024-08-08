KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

KOKKOSKERNELS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test_common)

if (KOKKOSKERNELS_ENABLE_EXPERIMENTAL)
  KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
          wiki_bsrmatrix
          SOURCES KokkosSparse_wiki_bsrmatrix.cpp
  )
endif()

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
        wiki_bsrmatrix_2
        SOURCES KokkosSparse_wiki_bsrmatrix_2.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  wiki_crsmatrix
  SOURCES KokkosSparse_wiki_crsmatrix.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  wiki_spmv
  SOURCES KokkosSparse_wiki_spmv.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  wiki_spadd
  SOURCES KokkosSparse_wiki_spadd.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  wiki_spgemm
  SOURCES KokkosSparse_wiki_spgemm.cpp
  )

KOKKOSKERNELS_ADD_EXECUTABLE_AND_TEST(
  wiki_gauss_seidel
  SOURCES KokkosSparse_wiki_gauss_seidel.cpp
  )
