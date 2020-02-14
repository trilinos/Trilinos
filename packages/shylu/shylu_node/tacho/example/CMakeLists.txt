
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(SOURCES "")

FILE(GLOB SOURCES *.cpp)

SET(LIBRARIES shyluchol)

IF (TACHO_HAVE_KOKKOS_TASK)
  #
  # Supernodes
  # 
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleCholSupernodesDeprecatedTaskScheduler
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleCholSupernodes.cpp
  #   COMM serial mpi
  # )
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleCholSupernodesDeprecatedTaskSchedulerMultiple
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleCholSupernodes.cpp
  #   COMM serial mpi
  # )
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

  #
  # Solver
  #
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleSolverDoubleDeprecatedTaskScheduler
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleSolver_double.cpp
  #   COMM serial mpi
  # )
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleSolverDoubleDeprecatedTaskSchedulerMultiple
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleSolver_double.cpp
  #   COMM serial mpi
  # )
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
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleSolverDoubleComplexDeprecatedTaskScheduler
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleSolver_dcomplex.cpp
  #   COMM serial mpi
  # )
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleSolverDoubleComplexDeprecatedTaskSchedulerMultiple
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleSolver_dcomplex.cpp
  #   COMM serial mpi
  # )
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
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleDenseByBlocksDeprecatedTaskScheduler
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleDenseByBlocks.cpp
  #   COMM serial mpi
  # )
  # TRIBITS_ADD_EXECUTABLE(
  #   Tacho_ExampleDenseByBlocksDeprecatedTaskSchedulerMultiple
  #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
  #   NOEXEPREFIX
  #   SOURCES Tacho_ExampleDenseByBlocks.cpp
  #   COMM serial mpi
  # )
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
    # TRIBITS_ADD_EXECUTABLE(
    #   Tacho_ExamplePerfTestDeprecatedTaskScheduler
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER
    #   NOEXEPREFIX
    #   SOURCES Tacho_ExamplePerfTest.cpp
    #   COMM serial mpi
    # )
    # TRIBITS_ADD_EXECUTABLE(
    #   Tacho_ExamplePerfTestDeprecatedTaskSchedulerMultiple
    #   TARGET_DEFINES -DTACHO_USE_DEPRECATED_TASKSCHEDULER_MULTIPLE
    #   NOEXEPREFIX
    #   SOURCES Tacho_ExamplePerfTest.cpp
    #   COMM serial mpi
    # )
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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUTacho_SimpleSparseTest_File
          SOURCE_FILES test.mtx
          SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
          DEST_DIR ${CMAKE_CURRENT_BINARY_DIR}
  )
ENDIF()
