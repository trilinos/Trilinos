
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluichol)

# CrsMatrix data structure examples
# =================================

TRIBITS_ADD_EXECUTABLE(
  crs_matrix_base
  SOURCES crs_matrix_base.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  crs_matrix_view
  SOURCES crs_matrix_view.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  crs_row_view
  SOURCES crs_row_view.cpp
  COMM serial mpi
  )

# DenseMatrix data structure examples
# ===================================

TRIBITS_ADD_EXECUTABLE(
  dense_matrix_base
  SOURCES dense_matrix_base.cpp
  COMM serial mpi
  )

# TPL library interface examples
# ==============================

TRIBITS_ADD_EXECUTABLE(
  scotch
  SOURCES scotch.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  qthread
  SOURCES qthread.cpp
  COMM serial mpi
  )

# Hier/Flat matrix conversion vice versa
# ======================================

TRIBITS_ADD_EXECUTABLE(
  crs_matrix_helper
  SOURCES crs_matrix_helper.cpp
  COMM serial mpi
  )

# Matrix partitioning
# ===================

TRIBITS_ADD_EXECUTABLE(
  partition
  SOURCES partition.cpp
  COMM serial mpi
  )

# Matrix base with non-built-in datatypes
# =======================================

TRIBITS_ADD_EXECUTABLE(
  crs_matrix_base_future
  SOURCES crs_matrix_base_future.cpp
  COMM serial mpi
  )

# Sequential for loop for debugging
# =================================

TRIBITS_ADD_EXECUTABLE(
  sequential_for
  SOURCES sequential_for.cpp
  COMM serial mpi
  )

# Execution Factory 
# =================

TRIBITS_ADD_EXECUTABLE(
  team_factory
  SOURCES team_factory.cpp
  COMM serial mpi
  )

# Matrix + Execution Policy
# =========================

TRIBITS_ADD_EXECUTABLE(
  crs_team_view
  SOURCES crs_team_view.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  crs_task_view
  SOURCES crs_task_view.cpp
  COMM serial mpi
  )

# Sparse BLAS level 1
# ===================

TRIBITS_ADD_EXECUTABLE(
  dot
  SOURCES dot.cpp
  COMM serial mpi
  )

# Sparse IChol 
# ============

TRIBITS_ADD_EXECUTABLE(
  ichol_right_unblocked
  SOURCES ichol_right_unblocked.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ichol_right_blocked
  SOURCES ichol_right_blocked.cpp
  COMM serial mpi
  )

TRIBITS_ADD_EXECUTABLE(
  ichol_right_by_blocks
  SOURCES ichol_right_by_blocks.cpp
  COMM serial mpi
  )

# Benchmark
# ==========

TRIBITS_ADD_EXECUTABLE(
 ichol_performance
 SOURCES ichol_performance.cpp
 COMM serial mpi
 )

TRIBITS_ADD_EXECUTABLE(
 ichol_graphviz
 SOURCES ichol_graphviz.cpp
 COMM serial mpi
 )


TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUICHOL_SimpleSparseTest_File
        SOURCE_FILES test.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
