
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#SET(GTEST_SOURCE_DIR /home/mndevec/work/Trilinos/packages/kokkos/tpls/gtest)

SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}../kokkos/tpls/gtest)

INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})

SET(GTEST_SOURCE_DIR ${${PARENT_PACKAGE_NAME}_SOURCE_DIR}/../kokkos/tpls/gtest)
INCLUDE_DIRECTORIES(${GTEST_SOURCE_DIR})
TRIBITS_ADD_LIBRARY(
  kokkoskernels_gtest
  HEADERS ${GTEST_SOURCE_DIR}/gtest/gtest.h
  SOURCES ${GTEST_SOURCE_DIR}/gtest/gtest-all.cc
  TESTONLY
  )
  
 TRIBITS_ADD_EXECUTABLE_AND_TEST(
  CrsMatrixUnitTest
  SOURCES TestCrsMatrix.cpp UnitTestMain.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "  FAILED  "
  TESTONLYLIBS  kokkoskernels_gtest
  )
  
  
# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  copyIntegers
  SOURCES copyIntegers.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  TESTONLYLIBS  kokkoskernels_gtest   
  )
  
# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  findRelOffset
  SOURCES findRelOffset.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  TESTONLYLIBS  kokkoskernels_gtest
  )
  
  # This test needs to be able to control when Kokkos gets initialized,
	# so it can't use the Teuchos unit test "main" .cpp file.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  MKL_sparseMatrixHandle
  SOURCES MKL_sparseMatrixHandle.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  TESTONLYLIBS  kokkoskernels_gtest  
  )

# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  replaceSumInto
  SOURCES replaceSumInto.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  TESTONLYLIBS  kokkoskernels_gtest
  )
  
# This test needs to be able to control when Kokkos gets initialized,
# so it can't use the Teuchos unit test "main" .cpp file.
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  replaceSumIntoLonger
  SOURCES replaceSumIntoLonger.cpp
  COMM serial mpi
  NUM_MPI_PROCS 1
  TESTONLYLIBS  kokkoskernels_gtest  
  )
    

#
# Disabled test temporarily; see #978.
# 
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   spgemm_unit_test
#   SOURCES UnitTestMain.cpp KokkosKernels_Test_SPGEMM.cpp
#   COMM serial mpi
#   ARGS "--mtx-file wathen_30_30.mtx"
#   NUM_MPI_PROCS 1
#   FAIL_REGULAR_EXPRESSION "  FAILED  "
#   TESTONLYLIBS  kokkoskernels_gtest
#   )

# TRIBITS_ADD_TEST(
#     spgemm_unit_test
#     NAME spgemm_sherman
#     NUM_MPI_PROCS 1
#     COMM serial mpi
# 	ARGS "--mtx-file sherman1.mtx"
# 	FAIL_REGULAR_EXPRESSION "  FAILED  "
# )

# TRIBITS_ADD_TEST(
#     spgemm_unit_test
#     NAME spgemm_si2
#     NUM_MPI_PROCS 1
#     COMM serial mpi
# 	ARGS "--mtx-file Si2.mtx"
# 	FAIL_REGULAR_EXPRESSION "  FAILED  "
# )

#
# Disabled temporarily per #979.
#

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   coloring_unit_test
#   SOURCES UnitTestMain.cpp KokkosKernels_Test_GraphColoring.cpp
#   COMM serial mpi
#   ARGS "--mtx-file wathen_30_30.mtx"
#   NUM_MPI_PROCS 1
#   FAIL_REGULAR_EXPRESSION "  FAILED  "
#   TESTONLYLIBS  kokkoskernels_gtest
#   )
  
# TRIBITS_ADD_TEST(
#     coloring_unit_test
#     NAME coloring_sherman
#     NUM_MPI_PROCS 1
#     COMM serial mpi
# 	ARGS "--mtx-file sherman1.mtx"
# 	FAIL_REGULAR_EXPRESSION "  FAILED  "
# )

# TRIBITS_ADD_TEST(
#     coloring_unit_test
#     NAME coloring_si2
#     NUM_MPI_PROCS 1
#     COMM serial mpi
# 	ARGS "--mtx-file Si2.mtx"
# 	FAIL_REGULAR_EXPRESSION "  FAILED  "
# )

#
# Disabled test temporarily; see #977.
# 
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   gauss_seidel_unit_test
#   SOURCES UnitTestMain.cpp KokkosKernels_Test_GaussSeidel.cpp
#   COMM serial mpi
#   ARGS "--mtx-file wathen_30_30.mtx"
#   NUM_MPI_PROCS 1
#   FAIL_REGULAR_EXPRESSION "  FAILED  "
#   TESTONLYLIBS  kokkoskernels_gtest
#   )

# TRIBITS_ADD_TEST(
#     gauss_seidel_unit_test
#     NAME gs_sherman
#     NUM_MPI_PROCS 1
#     COMM serial mpi
# 	ARGS "--mtx-file sherman1.mtx"
# 	FAIL_REGULAR_EXPRESSION "  FAILED  "
# )

#TRIBITS_ADD_TEST(
#    gauss_seidel_unit_test
#    NAME gs_si2
#    NUM_MPI_PROCS 1
#    COMM serial mpi
#	ARGS "--mtx-file Si2.mtx"
#	FAIL_REGULAR_EXPRESSION "  FAILED  "
#)
  
  TRIBITS_COPY_FILES_TO_BINARY_DIR(copy_files_for_kk_tests
  SOURCE_FILES
    sherman1.mtx
    Si2.mtx
    wathen_30_30.mtx
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
  DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
