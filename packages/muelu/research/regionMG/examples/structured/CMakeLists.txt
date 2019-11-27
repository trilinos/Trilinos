INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../../test/unit_tests)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../../mmayr/composite_to_regions/src)

# This test requires Tpetra and Amesos2, so it's only included if Tpetra is enabled
IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Amesos2)

  TRIBITS_ADD_EXECUTABLE(
    StructuredRegionDriver
    SOURCES Driver_Structured_Regions.cpp
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(StructuredRegion_cp
    SOURCE_FILES
      amg_1dof.xml
      structured_1dof.xml
      structured_1dof_3level.xml
      structured_1dof-complex.xml
      structured_linear_1dof.xml
      structured_linear_3dof.xml
      structured_unstructured_1dof.xml
      compare_residual_history.py
      gold_files/Star2D_GS_4.log.gold
      gold_files/Brick3D_linear_4.log.gold
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
    COMM serial mpi
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
    COMM serial mpi
    NUM_MPI_PROCS 4
    )
    
  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_CoarseSmoother_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --coarseSolverType=smoother"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_CoarseSmoother_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --coarseSolverType=smoother"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_Chebyshev"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_Chebyshev_1.log --smootherType=Chebyshev"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_Tpetra_Chebyshev"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_Chebyshev_4.log --smootherType=Chebyshev"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Star2D_AMG_CoarseSolver_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_1dof_3level.xml --matrixType=Star2D --nx=50 --ny=50 --smootherIts=2 --coarseSolverType=amg --coarseAmgXml=amg_1dof.xml"
    COMM serial mpi
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
    COMM serial mpi
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
    COMM serial mpi
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
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Unstructured_Laplace3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_unstructured_1dof.xml --matrixType=Laplace3D --nx=9 --ny=9 --nz=5 --smootherIts=2 --unstructured={3}"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=7 --ny=7 --nz=7 --smootherIts=2 --tol=1.0e-8 --convergence-log=Elasticity3D_linear_1.log"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

TRIBITS_ADD_TEST(
    StructuredRegionDriver
    NAME "Structured_Region_Linear_Elasticity3D_Tpetra"
    ARGS "--linAlgebra=Tpetra --xml=structured_linear_3dof.xml --matrixType=Elasticity3D --nx=7 --ny=7 --nz=7 --smootherIts=2 --tol=1.0e-8 --convergence-log=Elasticity3D_linear_4.log"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  ## These are very few regression tests, that also compare residual convergence history
  ASSERT_DEFINED(PYTHONINTERP_FOUND)
  IF (PYTHONINTERP_FOUND)

	TRIBITS_ADD_ADVANCED_TEST(
	  Structured_Region_Star2D_Tpetra_GS_MPI_4_regression
	  TEST_0
	    EXEC StructuredRegionDriver
 	    ARGS --linAlgebra=Tpetra --xml=structured_1dof.xml --matrixType=Star2D --nx=10 --ny=10 --smootherIts=2 --convergence-log=Star2D_GS_4.log --smootherType=Gauss
		NUM_MPI_PROCS 4
	  TEST_1
	    CMND ${PYTHON_EXECUTABLE}
	    ARGS compare_residual_history.py gold_files/Star2D_GS_4.log.gold Star2D_GS_4.log 1.0e-12
	    STANDARD_PASS_OUTPUT
	  FAIL_FAST
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
	    CMND ${PYTHON_EXECUTABLE}
	    ARGS compare_residual_history.py gold_files/Brick3D_linear_4.log.gold Brick3D_linear_4.log 1.0e-12
	    STANDARD_PASS_OUTPUT
	  FAIL_FAST
	  COMM serial mpi
	  OVERALL_NUM_MPI_PROCS 4
	  )

  ENDIF()

ENDIF()
