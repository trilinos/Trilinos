TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test/unit_tests)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)

# This test requires Tpetra and Amesos2, so it's only included if Tpetra is enabled
IF(${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    StructuredRegionDriver
    SOURCES Driver_Structured_Regions.cpp
    COMM serial mpi
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    StructuredRegion_cp
    SOURCE_FILES
      amg_1dof.xml
      amg_1dof_rebalance.xml
      amg_3dof_rebalance.xml
      coarse_smoother.xml
      compare_residual_history.py
      gold_files/Brick3D_linear_4.log.gold
      gold_files/Elasticity3D_CG_linear_4.log.gold
      gold_files/Elasticity3D_linear_4.log.gold
      gold_files/Elasticity3D_Richardson_linear_4.log.gold
      gold_files/Star2D_GS_4.log.gold
      gold_files/Star2D_W.log.gold
      structured_1dof-complex.xml
      structured_1dof.xml
      structured_1dof_3level.xml
      structured_linear_1dof.xml
      structured_linear_1dof_comp.xml
      structured_linear_3dof.xml
      structured_linear_R_1dof.xml
      structured_linear_R_3dof.xml
      structured_unstructured_1dof.xml
  )

  ## These are some completion tests
  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_GS"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_GS_1.log --smootherType=Gauss"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_GS"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_GS_4.log --smootherType=Gauss"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_SGS"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --smootherType=SymmetricGauss"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_SGS"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --smootherType=SymmetricGauss"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_CoarseSmoother_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --coarseSolverType=smoother --coarseSmootherXML=coarse_smoother.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_CoarseSmoother_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --coarseSolverType=smoother --coarseSmootherXML=coarse_smoother.xml"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_Chebyshev"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_Chebyshev_1.log --smootherType=Chebyshev --smootherChebyBoostFactor=1.2"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_Chebyshev"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_Chebyshev_4.log --smootherType=Chebyshev --smootherChebyBoostFactor=1.2"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --coarseSolverType=amg --coarseAmgXml=amg_1dof.xml"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_coarseRepartition_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --coarseSolverType=amg --coarseAmgXml=amg_1dof.xml --rebalance-coarse --keep-coarse-coords"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_repartitioning_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --coarseSolverType=amg --keep-coarse-coords --coarseAmgXml=amg_1dof_rebalance.xml"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_coarseRepartition_repartitioning_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --coarseSolverType=amg --keep-coarse-coords --coarseAmgXml=amg_1dof_rebalance.xml --rebalance-coarse"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Star2D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_linear_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_R_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_R_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_R_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_R_Brick3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_R_1dof.xml --matrixType=Brick3D --nx=19 --ny=19 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_R_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_R_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_R_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_linear_R_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_R_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_R_3dof.xml --matrixType=Elasticity3D --nx=19 --ny=19 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_linear_R_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Unstructured_Laplace3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_unstructured_1dof.xml --matrixType=Laplace3D --nx=9 --ny=9 --nz=5 --smootherIts=2 --unstructured={3}"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --tol=1.0e-8 --convergence-log=Elasticity3D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=7 --ny=7 --nz=7 --smootherIts=2 --tol=1.0e-8 --convergence-log=Elasticity3D_linear_4.log"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_repartitioning_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --tol=1.0e-8  --coarseSolverType=amg --keep-coarse-coords --coarseAmgXml=amg_3dof_rebalance.xml"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_coarseRepartition_repartitioning_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --tol=1.0e-8  --coarseSolverType=amg --keep-coarse-coords --rebalance-coarse --numPartitions=2 --coarseAmgXml=amg_3dof_rebalance.xml"
    COMM mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_Wcycle"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --cycleType=W"
    COMM serial mpi
    NUM_MPI_PROCS 1
  )

  ## These are regression tests, that also compare residual convergence history
  ASSERT_DEFINED(Python3_EXECUTABLE)
  IF(Python3_EXECUTABLE)

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Star2D_Tpetra_GS_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_GS_4.log --smootherType=Gauss
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Star2D_GS_4.log.gold Star2D_GS_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_Brick3D_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_1dof.xml --matrixType=Brick3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_4.log
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Brick3D_linear_4.log.gold Brick3D_linear_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_Elasticity3D_Region_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_linear_4.log --solverType=region
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Elasticity3D_linear_4.log.gold Elasticity3D_linear_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_Elasticity3D_CG_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_CG_linear_4.log --solverType=CG
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Elasticity3D_CG_linear_4.log.gold Elasticity3D_CG_linear_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_Elasticity3D_Richardson_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_Richardson_linear_4.log
           --solverType=Richardson
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Elasticity3D_Richardson_linear_4.log.gold Elasticity3D_Richardson_linear_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_R_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_1dof_comp.xml --matrixType=Brick3D --nx=19 --ny=19 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_4_comp.log
      NUM_MPI_PROCS 4
      TEST_1
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_R_1dof.xml --matrixType=Brick3D --nx=19 --ny=19 --nz=10 --smootherIts=2 --convergence-log=Brick3D_linear_R_4_comp.log
      NUM_MPI_PROCS 4
      TEST_2
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py Brick3D_linear_4_comp.log Brick3D_linear_R_4_comp.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Linear_Elasticity3D_Region_vs_Richardson_Tpetra_MPI_4_regression
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_Region_linear_4.log --solverType=region
      NUM_MPI_PROCS 4
      TEST_1
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=10 --ny=10 --nz=10 --smootherIts=2 --convergence-log=Elasticity3D_Richardson_linear_4.log
           --solverType=Richardson
      NUM_MPI_PROCS 4
      TEST_2
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py Elasticity3D_Region_linear_4.log Elasticity3D_Richardson_linear_4.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

    TRIBITS_ADD_ADVANCED_TEST(
      Structured_Region_Star2D_Tpetra_Wcycle
      TEST_0
      EXEC StructuredRegionDriver
      ARGS --linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --cycleType=W --convergence-log=Star2D_W.log
      NUM_MPI_PROCS 4
      TEST_1
      CMND ${Python3_EXECUTABLE}
      ARGS compare_residual_history.py gold_files/Star2D_W.log.gold Star2D_W.log 1.0e-12
      STANDARD_PASS_OUTPUT FAIL_FAST
      COMM serial mpi
      OVERALL_NUM_MPI_PROCS 4
    )

  ENDIF()

ENDIF()
