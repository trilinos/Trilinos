
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluchol)

# Dense Matrix Objects
# ====================

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixBase_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixBase_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixView_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixView_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_Cuda)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseMatrixBase_Cuda
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseMatrixBase_Cuda.cpp
    COMM serial mpi
  )
ENDIF()

# Dense Linear Algebra
# ====================

IF(Kokkos_ENABLE_OpenMP)
IF(TPL_ENABLE_MKL)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseGemmMKL_OpenMP
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseGemmMKL_OpenMP.cpp
    COMM serial mpi
  )
ENDIF()
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseGemmByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseGemmByBlocks_Pthreads.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseCholByBlocks_Pthreads
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseCholByBlocks_Pthreads.cpp
    COMM serial mpi
  )
ENDIF()


# # Do not compile if Scotch is not enabled
# IF(TPL_ENABLE_Scotch)

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_performance_device_cuda
#     NOEXEPREFIX
#     SOURCES example_chol_performance_device_cuda.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Trilinos_ENABLE_Teuchos)

# # symbolic factorization is guarded with Cholmod
# IF(TPL_ENABLE_Cholmod)

# # Symbolic Factorization
# # ======================

# IF(Kokkos_ENABLE_Serial)
#   TRIBITS_ADD_EXECUTABLE(
#     symbolic_factor_serial
#     NOEXEPREFIX
#     SOURCES example_symbolic_factor_serial.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     symbolic_factor_pthread
#     NOEXEPREFIX
#     SOURCES example_symbolic_factor_pthread.cpp
#     COMM serial mpi
#   )
# ENDIF()

# ENDIF()

# # Sparse Chol
# # ============

# IF(Kokkos_ENABLE_Serial)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_unblocked_serial
#     NOEXEPREFIX
#     SOURCES example_chol_unblocked_serial.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_by_blocks_graphviz
#     NOEXEPREFIX
#     SOURCES example_chol_by_blocks_graphviz.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_by_blocks_serial
#     NOEXEPREFIX
#     SOURCES example_chol_by_blocks_serial.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_unblocked_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_unblocked_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_by_blocks_pthread.cpp
#     COMM serial mpi
#   )
# ENDIF()

# # Triangular solve
# # ================

# IF(Kokkos_ENABLE_Serial)
#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_unblocked_serial
#     NOEXEPREFIX
#     SOURCES example_tri_solve_unblocked_serial.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_by_blocks_graphviz
#     NOEXEPREFIX
#     SOURCES example_tri_solve_by_blocks_graphviz.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_by_blocks_serial
#     NOEXEPREFIX
#     SOURCES example_tri_solve_by_blocks_serial.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_unblocked_pthread
#     NOEXEPREFIX
#     SOURCES example_tri_solve_unblocked_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_tri_solve_by_blocks_pthread.cpp
#     COMM serial mpi
#   )
# ENDIF()

# # Benchmark
# # ==========

# IF(Kokkos_ENABLE_Serial)
#   TRIBITS_ADD_EXECUTABLE(
#     stat_by_blocks_serial
#     NOEXEPREFIX
#     SOURCES example_stat_by_blocks_serial.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_performance_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_performance_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_performance_single_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_performance_single_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_performance_pthread
#     NOEXEPREFIX
#     SOURCES example_tri_solve_performance_pthread.cpp
#     COMM serial mpi
#   )

#   # TRIBITS_ADD_EXECUTABLE(
#   #   kokkos_task_data_pthread
#   #   NOEXEPREFIX
#   #   SOURCES example_kokkos_task_data_pthread.cpp
#   #   COMM serial mpi
#   # )

#   # TRIBITS_ADD_EXECUTABLE(
#   #   kokkos_data_data_pthread
#   #   NOEXEPREFIX
#   #   SOURCES example_kokkos_data_data_pthread.cpp
#   #   COMM serial mpi
#   # )

#   TRIBITS_ADD_EXECUTABLE(
#     kokkos_team_barrier_pthread
#     NOEXEPREFIX
#     SOURCES example_kokkos_team_barrier_pthread.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_QTHREAD)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_performance_qthread
#     NOEXEPREFIX
#     SOURCES example_chol_performance_qthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     tri_solve_performance_qthread
#     NOEXEPREFIX
#     SOURCES example_tri_solve_performance_qthread.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(TPL_ENABLE_HYPRE)
#   TRIBITS_ADD_EXECUTABLE(
#     hypre_euclid_performance
#     NOEXEPREFIX
#     SOURCES example_hypre_euclid_performance.cpp
#     COMM serial mpi
#   )
# ENDIF()

# # Pardiso
# # =======
# IF(Kokkos_ENABLE_OpenMP)
# IF(TPL_ENABLE_MKL)
#   TRIBITS_ADD_EXECUTABLE(
#     pardiso_chol_performance_openmp
#     NOEXEPREFIX
#     SOURCES example_pardiso_chol_performance_openmp.cpp
#     COMM serial mpi
#   )
# ENDIF()
# ENDIF()

# # Complete Cholesky
# # =================
# IF(Kokkos_ENABLE_Serial)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_direct_plain_serial
#     NOEXEPREFIX
#     SOURCES example_chol_direct_plain_serial.cpp
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     chol_direct_plain_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_direct_plain_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_direct_nested_dense_block_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_direct_nested_dense_block_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     chol_direct_nested_dense_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_chol_direct_nested_dense_by_blocks_pthread.cpp
#     COMM serial mpi
#   )

#   # TRIBITS_ADD_EXECUTABLE(
#   #   chol_direct_solver_pthread
#   #   NOEXEPREFIX
#   #   SOURCES example_chol_direct_solver_pthread.cpp
#   #   COMM serial mpi
#   # )
# ENDIF()

# # Dense Linear Algebra
# # ====================
# IF(Kokkos_ENABLE_Pthread)
#   TRIBITS_ADD_EXECUTABLE(
#     dense_gemm_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_dense_gemm_by_blocks_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     dense_trsm_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_dense_trsm_by_blocks_pthread.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     dense_chol_by_blocks_pthread
#     NOEXEPREFIX
#     SOURCES example_dense_chol_by_blocks_pthread.cpp 
#     COMM serial mpi
#   )
# ENDIF()

# IF(Kokkos_ENABLE_OpenMP)
# IF(TPL_ENABLE_MKL)
#   TRIBITS_ADD_EXECUTABLE(
#     dense_gemm_mkl_openmp
#     NOEXEPREFIX
#     SOURCES example_dense_gemm_mkl_openmp.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     dense_trsm_mkl_openmp
#     NOEXEPREFIX
#     SOURCES example_dense_trsm_mkl_openmp.cpp
#     COMM serial mpi
#   )

#   TRIBITS_ADD_EXECUTABLE(
#     dense_chol_mkl_openmp
#     NOEXEPREFIX
#     SOURCES example_dense_chol_mkl_openmp.cpp
#     COMM serial mpi
#   )
# ENDIF()

# ENDIF()  # Teuchos

# ENDIF()  # Scotch

# # Temporary
# # =========

# # TRIBITS_ADD_EXECUTABLE(
# #  kokkos_task_data
# #  NOEXEPREFIX
# #  SOURCES kokkos_task_data.cpp
# #  COMM serial mpi
# #  )

# ENDIF()

# TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUTacho_SimpleSparseTest_File
#         SOURCE_FILES test.mtx
#         SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
#         DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
# )
