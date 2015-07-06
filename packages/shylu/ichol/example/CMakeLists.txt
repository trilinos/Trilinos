
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluichol)

# Symbolic Factorization
# ======================

IF(Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE(
    symbolic_factor_serial
    SOURCES example_symbolic_factor_serial.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    symbolic_factor_pthread
    SOURCES example_symbolic_factor_pthread.cpp
    COMM serial mpi
  )
ENDIF()

# Sparse IChol
# ============

IF(Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE(
    ichol_unblocked_serial
    SOURCES example_ichol_unblocked_serial.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ichol_by_blocks_graphviz
    SOURCES example_ichol_by_blocks_graphviz.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ichol_by_blocks_serial
    SOURCES example_ichol_by_blocks_serial.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    ichol_unblocked_pthread
    SOURCES example_ichol_unblocked_pthread.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ichol_by_blocks_pthread
    SOURCES example_ichol_by_blocks_pthread.cpp
    COMM serial mpi
  )
ENDIF()

# Triangular solve
# ================

IF(Kokkos_ENABLE_Serial)
  TRIBITS_ADD_EXECUTABLE(
    tri_solve_unblocked_serial
    SOURCES example_tri_solve_unblocked_serial.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    tri_solve_by_blocks_graphviz
    SOURCES example_tri_solve_by_blocks_graphviz.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    tri_solve_by_blocks_serial
    SOURCES example_tri_solve_by_blocks_serial.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    tri_solve_unblocked_pthread
    SOURCES example_tri_solve_unblocked_pthread.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    tri_solve_by_blocks_pthread
    SOURCES example_tri_solve_by_blocks_pthread.cpp
    COMM serial mpi
  )
ENDIF()

# Benchmark
# ==========

IF(Kokkos_ENABLE_Pthread)
  TRIBITS_ADD_EXECUTABLE(
    ichol_performance_pthread
    SOURCES example_ichol_performance_pthread.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    tri_solve_performance_pthread
    SOURCES example_tri_solve_performance_pthread.cpp
    COMM serial mpi
  )
ENDIF()

IF(Kokkos_ENABLE_QTHREAD)
  TRIBITS_ADD_EXECUTABLE(
    ichol_performance_qthread
    SOURCES example_ichol_performance_qthread.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    tri_solve_performance_qthread
    SOURCES example_tri_solve_performance_qthread.cpp
    COMM serial mpi
  )
ENDIF()

# Temporary
# =========

# TRIBITS_ADD_EXECUTABLE(
#  kokkos_task_data
#  SOURCES kokkos_task_data.cpp
#  COMM serial mpi
#  )



TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUICHOL_SimpleSparseTest_File
        SOURCE_FILES test.mtx
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
        DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
)
