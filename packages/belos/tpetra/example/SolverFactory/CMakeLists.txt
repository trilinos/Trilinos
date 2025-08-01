ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Galeri)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Xpetra)

IF (${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Xpetra)

  TRIBITS_ADD_EXECUTABLE(
    SolverFactory_Tpetra_Galeri_Ex
    SOURCES SolverFactoryTpetraGaleriEx.cpp
    COMM serial mpi
    CATEGORIES BASIC PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    SolverFactory_Tpetra_Galeri_Ex
    NAME SolverFactory_Tpetra_Galeri_Ex
    COMM serial mpi
    ARGS --verbose
    STANDARD_PASS_OUTPUT
  )

  set(COUNTER 1)
  set(MAX_COUNT 11)

  WHILE(COUNTER LESS MAX_COUNT)
    math(EXPR NP "${COUNTER} * ${COUNTER}")

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_PseudoBlockCG
      COMM serial mpi
      ARGS "--verbose --solverName=\"Pseudoblock CG\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_BlockCG
      COMM serial mpi
      ARGS "--verbose --solverName=\"Block CG\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_SingleReduceCG
      COMM serial mpi
      ARGS "--verbose --solverName=\"Single-Reduce CG\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_PseudoBlockGMRES
      COMM serial mpi
      ARGS "--verbose --solverName=\"Pseudoblock GMRES\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_BlockGMRES
      COMM serial mpi
      ARGS "--verbose --solverName=\"Block GMRES\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    TRIBITS_ADD_TEST(
      SolverFactory_Tpetra_Galeri_Ex
      NAME PerfTest_BiCGStab
      COMM serial mpi
      ARGS "--verbose --solverName=\"BiCGStab\" --matrixType=Laplace3D --nx=200 --ny=200 --nz=200 --max-iters=100"
      NUM_MPI_PROCS ${NP}
      STANDARD_PASS_OUTPUT
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )

    math(EXPR COUNTER "${COUNTER} + 1")
  ENDWHILE()


ENDIF (${PACKAGE_NAME}_ENABLE_Galeri AND ${PACKAGE_NAME}_ENABLE_Xpetra)
