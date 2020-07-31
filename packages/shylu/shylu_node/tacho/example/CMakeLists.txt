
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES tacho)

TRIBITS_ADD_EXECUTABLE(
  Tacho_ExampleCombineDataFileToMatrixMarketFile
  NOEXEPREFIX
  SOURCES Tacho_ExampleCombineDataFileToMatrixMarketFile.cpp
  COMM serial mpi
)

IF (TACHO_HAVE_KOKKOS_TASK)
  #
  # Supernodes
  # 
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholSupernodesTaskScheduler
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholSupernodes.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholSupernodesTaskSchedulerMultiple
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholSupernodes.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholSupernodesChaseLevTaskScheduler
    TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholSupernodes.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholTriSolveVar0
    TARGET_DEFINES -DTACHO_USE_TRISOLVE_VARIANT=0
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholTriSolve.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholTriSolveVar1
    TARGET_DEFINES -DTACHO_USE_TRISOLVE_VARIANT=1
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholTriSolve.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholLevelSetVar0
    TARGET_DEFINES -DTACHO_USE_TRISOLVE_VARIANT=0
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholLevelSet.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleCholLevelSetVar1
    TARGET_DEFINES -DTACHO_USE_TRISOLVE_VARIANT=1
    NOEXEPREFIX
    SOURCES Tacho_ExampleCholLevelSet.cpp
    COMM serial mpi
  )

  #
  # Solver
  #
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleTaskScheduler
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_double.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleTaskSchedulerMultiple
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_double.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleChaseLevTaskScheduler
    TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_double.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleComplexTaskScheduler
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_dcomplex.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleComplexTaskSchedulerMultiple
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_dcomplex.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleSolverDoubleComplexChaseLevTaskScheduler
    TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleSolver_dcomplex.cpp
    COMM serial mpi
  )

  #
  # External Interface
  # 
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleExternalInterface
    NOEXEPREFIX
    SOURCES Tacho_ExampleExternalInterface.cpp
    COMM serial mpi
  )

  #
  # DenseByBlocks 
  #
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseByBlocksTaskScheduler
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseByBlocks.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseByBlocksTaskSchedulerMultiple
    TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseByBlocks.cpp
    COMM serial mpi
  )
  TRIBITS_ADD_EXECUTABLE(
    Tacho_ExampleDenseByBlocksChaseLevTaskScheduler
    TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
    NOEXEPREFIX
    SOURCES Tacho_ExampleDenseByBlocks.cpp
    COMM serial mpi
  )

  #
  # NVIDIA cuSolver
  # 
  IF(Kokkos_ENABLE_CUDA) 
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExampleCholCuSolver
      NOEXEPREFIX
      SOURCES Tacho_ExampleCholCuSolver.cpp
      COMM serial mpi
    )
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExampleCuSparseTriSolve
      NOEXEPREFIX
      SOURCES Tacho_ExampleCuSparseTriSolve.cpp
      COMM serial mpi
    )
  ENDIF()



  #
  # Intel MKL Pardiso and PerfTest
  #
  IF(Kokkos_ENABLE_OPENMP OR Kokkos_ENABLE_SERIAL)
    IF(TPL_ENABLE_MKL)
      TRIBITS_ADD_EXECUTABLE(
        Tacho_ExampleCholPardiso
        NOEXEPREFIX
        SOURCES Tacho_ExampleCholPardiso.cpp
        COMM serial mpi
      )
    ENDIF()
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExamplePerfTestTaskScheduler
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_ExamplePerfTest.cpp
      COMM serial mpi
    )
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExamplePerfTestTaskSchedulerMultiple
      TARGET_DEFINES -DTACHO_USE_TASKSCHEDULER_MULTIPLE
      NOEXEPREFIX
      SOURCES Tacho_ExamplePerfTest.cpp
      COMM serial mpi
    )
    TRIBITS_ADD_EXECUTABLE(
      Tacho_ExamplePerfTestChaseLevTaskScheduler
      TARGET_DEFINES -DTACHO_USE_CHASELEV_TASKSCHEDULER
      NOEXEPREFIX
      SOURCES Tacho_ExamplePerfTest.cpp
      COMM serial mpi
    )
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Tacho_SimpleSparseTest_File
          SOURCE_FILES test.mtx graph.dat weight.dat graph_test.mtx 
          SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
          DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()
