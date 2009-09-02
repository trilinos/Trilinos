
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(PackageCopyFilesToBinaryDir)
INCLUDE(MultilineSet)


ASSERT_DEFINED(PACKAGE_SOURCE_DIR CMAKE_CURRENT_SOURCE_DIR)


FILE(GLOB SOURCES *.cpp)


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NOX)
IF (${PACKAGE_NAME}_ENABLE_NOX)
  ADD_DEFINITIONS(-DHAVE_RYTHMOS_NOX)
ENDIF()

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Stratimikos)
IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
  ADD_DEFINITIONS(-DHAVE_RYTHMOS_STRATIMIKOS)
ENDIF()

PACKAGE_ADD_EXECUTABLE(
  complicatedExample
  SOURCES ${SOURCES}
  COMM serial mpi
  )


PACKAGE_ADD_TEST(
  complicatedExample
  NAME complicatedExample_explicit
  ARGS
    # Default run, Explicit Runge-Kutta 4 stage
    "--verbose --echo-command-line"
    # Scalar Forward Euler
    "--method=FE --maxerror=1e-1 --verbose --echo-command-line"
    # Vector Explicit Runge-Kutta 4 stage
    "--numelements=10 --verbose --echo-command-line"
    # Vector Forward Euler
    "--method=FE --maxerror=1e-1 --numelements=10 --verbose --echo-command-line"
    # Vector Forward Euler with same accuracy as ERK
    "--numelements=100 --method=FE --T=0.1 --numsteps=5000 --verbose --echo-command-line"
  COMM serial mpi  
  NUM_MPI_PROCS 1
  )


ASSERT_DEFINED(Stratimikos_ENABLE_Amesos)
IF (Stratimikos_ENABLE_Amesos)

  MULTILINE_SET( AMESOS_PL
    "<ParameterList>"
    "  <Parameter name='Linear Solver Type' type='string' value='Amesos'/>"
    "</ParameterList>"
    )

  PACKAGE_ADD_TEST(
    complicatedExample
    NAME complicatedExample_amesos
    ARGS
      # Scalar Backward Euler
      "--method=BE --maxerror=1e-1 --verbose --echo-command-line"
      "--method=BDF --numelements=1 --x0=1 --lambda_fit=linear --lambda_min=-0.3 --lambda_max=-0.3 --stepmethod=variable --T=1 --force_coeff=31.4159 --maxorder=1 --reltol=1e-3 --abstol=1e-3 --maxerror=1e-3 --max-restep-error=1e-9  --verbose --echo-command-line"
      # Vector Backward Euler
      "--method=BE --maxerror=1e-1 --numelements=10 --verbose --echo-command-line"
      # Scalar Backward Euler using Amesos
      "--method=BE --maxerror=1e-1 --numelements=1 --verbose --extra-linear-solver-params=\"${AMESOS_PL}\" --linear-solver-params-used-file=lowsf.amesos.used.xml --echo-command-line"
      # Vector Backward Euler using Amesos
      "--method=BE --maxerror=1e-1 --numelements=10 --verbose --extra-linear-solver-params=\"${AMESOS_PL}\" --linear-solver-params-used-file=lowsf.amesos.used.xml --echo-command-line --verb-level=VERB_MEDIUM"
      # Vector Backward Euler with lower accuracy than ERK
      "--numelements=10 --method=BE --T=0.1 --numsteps=500 --maxerror=1e-5 --verbose --echo-command-line --verb-level=VERB_MEDIUM"
      # Vector Implicit BDF
      "--numelements=10 --method=BDF --stepmethod=variable --reltol=1e-6 --abstol=1e-6 --T=0.1 --maxerror=1e-4 --verbose --echo-command-line"
    COMM serial mpi  
    NUM_MPI_PROCS 1
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NOX)
  IF (${PACKAGE_NAME}_ENABLE_NOX)
    PACKAGE_ADD_TEST(
      complicatedExample
      NAME complicatedExample_amesos
      POSTFIX_AND_ARGS_0 nox_bdf
        --numelements=10 --method=BDF --stepmethod=variable --reltol=1e-6
        --abstol=1e-6 --T=0.1 --maxerror=1e-4 --useNOX --verbose --echo-command-line
      POSTFIX_AND_ARGS_1 nox_bd
        --method=BDF --maxorder=1 --maxerror=1e-1 --numelements=10 --verbose
        --echo-command-line
      COMM serial mpi  
      NUM_MPI_PROCS 1
      )
  ENDIF()


ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_AztecOO)
IF (Stratimikos_ENABLE_AztecOO)

  MULTILINE_SET(AZTECOO_PL
    "<ParameterList>"
    "  <Parameter name='Linear Solver Type' type='string' value='AztecOO'/>"
    "  <ParameterList name ='Linear Solver Types'>"
    "    <ParameterList name='AztecOO'>"
    "      <ParameterList name='Forward Solve'>"
    "        <Parameter name='Tolerance' type='double' value='1e-12'/>"
    "      </ParameterList>"
    "    </ParameterList>"
    "  </ParameterList>"
    "</ParameterList>"
    )

  PACKAGE_ADD_TEST(
    complicatedExample
    NAME complicatedExample_aztecoo
    ARGS
      # Scalar Backward Euler using AztecOO
      "--method=BE --maxerror=1e-1 --numelements=1 --verbose --extra-linear-solver-params=\"${AZTECOO_PL}\" --linear-solver-params-used-file=lowsf.aztecoo.used.xml --echo-command-line"
      # Vector Backward Euler using AztecOO
      "--method=BE --maxerror=1e-1 --numelements=100 --verbose --extra-linear-solver-params=\"${AZTECOO_PL}\" --linear-solver-params-used-file=lowsf.aztecoo.used.xml --echo-command-line"
      # Vector Implicit BDF with AztecOO
      "--numelements=100 --method=BDF --stepmethod=variable --reltol=1e-6 --abstol=1e-6 --T=0.1 --maxerror=1e-4 --verbose --extra-linear-solver-params=\"${AZTECOO_PL}\" --linear-solver-params-used-file=lowsf.aztecoo.used.xml --echo-command-line"
    COMM serial mpi  
    NUM_MPI_PROCS 1
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NOX)
  IF (${PACKAGE_NAME}_ENABLE_NOX)

    PACKAGE_ADD_TEST(
      complicatedExample
      NAME complicatedExample_aztecoo
      POSTFIX_AND_ARGS_0 nox
        --numelements=100 --method=BDF --stepmethod=variable --reltol=1e-6 --abstol=1e-6
        --T=0.1 --maxerror=1e-4 --useNOX --verbose --extra-linear-solver-params="${AZTECOO_PL}"
        --linear-solver-params-used-file=lowsf.aztecoo.used.xml --echo-command-line
      COMM serial mpi  
      NUM_MPI_PROCS 1
      )

  ENDIF()

ENDIF()


ASSERT_DEFINED(Stratimikos_ENABLE_Belos)
IF (Stratimikos_ENABLE_Belos)

  MULTILINE_SET(BELOS_BLOCK_GMRES_PL
    "<ParameterList>"
    "  <Parameter name='Linear Solver Type' type='string' value='Belos'/>"
    "  <ParameterList name='Linear Solver Types'>"
    "    <ParameterList name='Belos'>"
    "      <Parameter name='Solver Type' type='string' value='Block GMRES'/>"
    "      <ParameterList name='Solver Types'>"
    "        <ParameterList name='Block GMRES'>"
    "          <Parameter name='Convergence Tolerance' type='double' value='1e-12'/>"
    "        </ParameterList>"
    "      </ParameterList>"
    "    </ParameterList>"
    "  </ParameterList>"
    "</ParameterList>"
    )


  PACKAGE_ADD_TEST(
    complicatedExample
    NAME complicatedExample_belos
    ARGS
      # Vector Backward Euler using Belos
      "--method=BE --maxerror=1e-1 --numelements=100 --verbose --extra-linear-solver-params=\"${BELOS_BLOCK_GMRES_PL}\" --linear-solver-params-used-file=lowsf.belos.used.xml"
      # Vector Implicit BDF with Belos
      "--numelements=100 --method=BDF --stepmethod=variable --reltol=1e-6 --abstol=1e-6 --T=0.1 --maxerror=1e-4 --verbose --extra-linear-solver-params=\"${BELOS_BLOCK_GMRES_PL}\" --linear-solver-params-used-file=lowsf.belos.used.xml" 
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_NOX)
  IF (${PACKAGE_NAME}_ENABLE_NOX)
    PACKAGE_ADD_TEST(
      complicatedExample
      NAME complicatedExample_belos
      POSTFIX_AND_ARGS_0 nox_bdf
        --numelements=100 --method=BDF --stepmethod=variable --reltol=1e-6 
        --abstol=1e-6 --T=0.1 --maxerror=1e-4 --useNOX --verbose
        --extra-linear-solver-params="${BELOS_BLOCK_GMRES_PL}"
        --linear-solver-params-used-file=lowsf.belos.used.xml 
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  ENDIF()

ENDIF()
