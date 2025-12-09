# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Amesos
  ${PACKAGE_NAME}_ENABLE_Amesos2
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_Epetra
)

#
# Executable
#

#TODO: how to remove that?
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

IF(${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  SET(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK YES)
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra
   AND ${PACKAGE_NAME}_ENABLE_EpetraExt
   AND ${PACKAGE_NAME}_ENABLE_Ifpack
   AND ${PACKAGE_NAME}_ENABLE_Amesos
)
  SET(${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK YES)
ENDIF()

IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK OR ${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK)

  TRIBITS_ADD_EXECUTABLE(
    Driver
    SOURCES Driver.cpp
    COMM serial mpi
    CATEGORIES BASIC PERFORMANCE
  )

  INSTALL(TARGETS "${PACKAGE_NAME}_Driver")

  TRIBITS_ADD_TEST(
    Driver
    NAME Defaults
    COMM serial mpi
    ARGS "--xml=DEFAULTS"
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  # Perf Model
  TRIBITS_ADD_TEST(
    Driver
    NAME PerformanceModel
    COMM mpi
    ARGS "--nx=40 --ny=40 --nz=40 --matrixType=Laplace3D --performance-model=verbose --its=14"
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  # Do a simple weak scaling experiment in 3D
  SET(COUNTER 1)
  SET(MAX_COUNT 7)
  WHILE(COUNTER LESS MAX_COUNT)
    MATH(EXPR NP "${COUNTER} * ${COUNTER} * ${COUNTER}")
    MATH(EXPR NX "${COUNTER} * 120")
    TRIBITS_ADD_TEST(
      Driver
      NAME SetupSolve_Performance
      COMM mpi
      ARGS "--stacked-timer --nx=${NX} --ny=${NX} --nz=${NX} --matrixType=Laplace3D"
      NUM_MPI_PROCS ${NP}
      PASS_REGULAR_EXPRESSION "Belos converged"
      RUN_SERIAL
      CATEGORIES PERFORMANCE
    )
    MATH(EXPR COUNTER "${COUNTER} + 1")
  ENDWHILE()

  # This an artifact to reproduce the old 4-rank test that is weak scaled by only adding guys in (x,y)
  TRIBITS_ADD_TEST(
    Driver
    NAME SetupSolve_Performance
    COMM mpi
    ARGS "--stacked-timer --nx=240 --ny=240 --nz=120 --matrixType=Laplace3D"
    NUM_MPI_PROCS 4
    PASS_REGULAR_EXPRESSION "Belos converged"
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    Driver_cp
    SOURCE_FILES
      amgx.json
      amgx.xml
      anisotropic.xml
      aux1.mat
      aux2.mat
      blkDiag3x3And6x6.mat
      blkSmooEquivOlapSchwarz.xml
      block_diagonal_drop.xml
      blocknumber.mm
      circ_nsp_dependency.xml
      comboP.xml
      comp_rotations.xml
      conchas_milestone_zoltan.xml
      conchas_milestone_zoltan2.xml
      conchas_milestone_zoltan2_complex.xml
      coords1.mat
      coords2.mat
      coords2.mm
      elasticity3D.xml
      emin_with_dl_elasticity3D.xml
      emin_with_dl_material_elasticity3D.xml
      GblkAmat.dat
      GblkMap.dat
      GblkRhs.dat
      Gblks.dat
      generalBlkSmoothing.xml
      iso_poisson.xml
      isorropia.xml
      material.mm
      material.xml
      multiPhys4x4.mat
      multiPhys4x4.mat
      null1.mat
      null2.mat
      oLapSchwarzEquivBlkSmoo.xml
      pairwise.xml
      rap.xml
      recurMG.xml
      regularOverLap.dat
      repart.xml
      sa_enforce_constraints.xml
      sa_with_Ifpack2_line_detection.xml
      sa_with_ilu.xml
      scaling-complex-withglobalconstants.xml
      scaling-complex.xml
      scaling-gpu.xml
      scaling-ml.xml
      scaling-with-rerun.xml
      scaling-withglobalconstants.xml
      scaling.xml
      scaling.yaml
      scaling_apex.xml
      scaling_distance2_agg.xml
      scaling_sendtype.xml
      smoother.xml
      smoother_complex.xml
      smooVec.mm
      smooVecCoalesce.xml
      tripleMatrixProduct.xml
      zeroAggTest.xml
    CATEGORIES BASIC PERFORMANCE
  )

  TRIBITS_ADD_EXECUTABLE(
    ImportPerformance
    SOURCES ImportPerformance.cpp
    COMM mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ComboPDriver
    SOURCES ComboPDriver.cpp
    COMM mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    TAFCPerformance
    SOURCES TAFCPerformance.cpp
    COMM mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    ReadMatrix
    SOURCES ReadMatrix.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_EXECUTABLE(
    MatrixMatrixMultiply
    SOURCES MatrixMatrixMultiply.cpp
    COMM mpi
    CATEGORIES BASIC PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    MatrixMatrixMultiply
    NAME MatrixMultiply_Performance_1
    COMM mpi
    ARGS "--timings --seed=12345 --minrows=8000 --maxrows=10000 --nmults=5"
    NUM_MPI_PROCS 1
    STANDARD_PASS_OUTPUT
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    MatrixMatrixMultiply
    NAME MatrixMultiply_Performance_4
    COMM mpi
    ARGS "--timings --seed=12345 --minrows=8000 --maxrows=10000 --nmults=5"
    NUM_MPI_PROCS 4
    STANDARD_PASS_OUTPUT
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME Driver_anisotropic
    #ARGS "--xml=anisotropic.xml --tol=1e-6 --its=66"
    #FIXME JHU: 2023-Nov-20 In OpenMP build, this test converges in 190 iterations
    #FIXME                  For other builds, test converges in less than 66 iterations
    ARGS "--xml=anisotropic.xml --tol=1e-6"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME Driver_sendtype
    ARGS "--xml=scaling_sendtype.xml"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  # not very elegant. The best would probably to make Reuse.cpp work for both Epetra and Tpetra
  IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)
    TRIBITS_ADD_EXECUTABLE(
      Reuse
      SOURCES Reuse.cpp
      COMM serial mpi
    )
    TRIBITS_COPY_FILES_TO_BINARY_DIR(Reuse_cp SOURCE_FILES reuse_emin_emin.xml reuse_sa_RP.xml)
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    ReadMatrix_cp
    SOURCE_FILES
      A.mm
      B.mm
      coords.mm
      simple.xml
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(ImportPerformance_cp SOURCE_FILES import.xml)

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    ComboP_cp
    SOURCE_FILES
      aux1.mat
      aux2.mat
      comboP.xml
      coords1.mat
      coords2.mat
      multiPhys4x4.mat
      multiPhys4x4.mat
      null1.mat
      null2.mat
  )

ENDIF()

IF(${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK)

  TRIBITS_ADD_EXECUTABLE(
    ImportPerformanceEpetra
    SOURCES ImportPerformance
    COMM mpi
    TARGET_DEFINES ${PACKAGE_NAME}_UNDEFINE_Tpetra=1
  )

ENDIF()

#
# Tests
#

IF(${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK AND (NOT Xpetra_INT_LONG_LONG))

  IF(${PACKAGE_NAME}_ENABLE_Zoltan)
    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_Defaults"
      ARGS "--linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
    )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_Rerun"
      ARGS "--linAlgebra=Epetra --xml=scaling-with-rerun.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
    )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_isotropic_poisson"
      ARGS "--linAlgebra=Epetra --xml=iso_poisson.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
    )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverEpetra_Milestone"
      ARGS "--linAlgebra=Epetra --xml=conchas_milestone_zoltan.xml"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
    )

    IF(${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_Isorropia)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverEpetra_IsorropiaPoisson"
        ARGS "--linAlgebra=Epetra --xml=isorropia.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ENDIF()

    IF(${PACKAGE_NAME}_ENABLE_Zoltan AND ${PACKAGE_NAME}_ENABLE_Isorropia)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverEpetra_CircNspDependency"
        ARGS "--linAlgebra=Epetra --xml=circ_nsp_dependency.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ENDIF()

  ENDIF()

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "RAPScalingTestEpetra"
    ARGS "--linAlgebra=Epetra --xml=rap.xml --solver=none --nx=50 --ny=50 --rebuild=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "SmootherScalingTestEpetra"
    ARGS "--linAlgebra=Epetra --xml=smoother.xml --nx=10 --ny=10 --solver=none"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    ReadMatrix
    NAME "ReadMatrixEpetra"
    ARGS "--linAlgebra=Epetra --xml=simple.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
  )

ENDIF()

IF(${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK)

  TRIBITS_ADD_TEST(
    ImportPerformanceEpetra
    NAME "ImportPerformance_Epetra"
    ARGS "--linAlgebra=Epetra --nx=30 --ny=30"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

ENDIF()

IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)

  TRIBITS_ADD_TEST(
    ImportPerformance
    NAME "ImportPerformance_Tpetra"
    ARGS "--linAlgebra=Tpetra --nx=30 --ny=30"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    #do not check for convergence
  )

  TRIBITS_ADD_TEST(
    ComboPDriver
    NAME "ComboPTest"
    ARGS "--linAlgebra=Tpetra --its=9"
    NUM_MPI_PROCS 4
    COMM mpi
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

ENDIF()

IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)

  TRIBITS_ADD_TEST(
    Driver
    NAME "GeneralBlockSmoothing"
    ARGS "--linAlgebra=Tpetra --xml=generalBlkSmoothing.xml --belosType=\"Fixed\ Point\" --rowmap=GblkMap.dat --matrix=GblkAmat.dat --rhs=GblkRhs.dat --its=2 --userBlks=Gblks.dat --notimings --tol=.976 --its=2"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )
  TRIBITS_ADD_TEST(
    Driver
    NAME "SerialGeneralBlockSmoothing"
    ARGS "--linAlgebra=Tpetra --xml=generalBlkSmoothing.xml --belosType=\"Fixed\ Point\" --matrix=GblkAmat.dat --rhs=GblkRhs.dat --its=2 --userBlks=Gblks.dat --notimings --tol=.976 --its=2"
    NUM_MPI_PROCS 1
    COMM mpi # HAVE_MPI required
    RUN_SERIAL
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "ZeroAggTest"
    ARGS "--linAlgebra=Tpetra --xml=zeroAggTest.xml --belosType=\"Block GMRES\" --matrix=blkDiag3x3And6x6.mat --its=11 --notimings --coords=coords2.mm"
    NUM_MPI_PROCS 2
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )
  TRIBITS_ADD_TEST(
    Driver
    NAME "SerialZeroAggTest"
    ARGS "--linAlgebra=Tpetra --xml=zeroAggTest.xml --belosType=\"Block GMRES\" --matrix=blkDiag3x3And6x6.mat --its=11 --notimings --coords=coords2.mm"
    NUM_MPI_PROCS 1
    COMM mpi # HAVE_MPI required
    RUN_SERIAL
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "BlockSmoothingWithAverages"
    ARGS "--linAlgebra=Tpetra --xml=blkSmooEquivOlapSchwarz.xml --belosType=\"Fixed\ Point\" --rowmap=GblkMap.dat --matrix=GblkAmat.dat --rhs=GblkRhs.dat --tol=.1 --userBlks=regularOverLap.dat --its=12"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "SchwarzSmoothingWithAverages"
    ARGS "--linAlgebra=Tpetra --xml=oLapSchwarzEquivBlkSmoo.xml --belosType=\"Fixed\ Point\" --rowmap=GblkMap.dat --matrix=GblkAmat.dat --rhs=GblkRhs.dat --tol=.1 --its=12"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "MaterialDistanceLaplacian"
    POSTFIX_AND_ARGS_0 "ConstantScalar" --linAlgebra=Tpetra --xml=material.xml --matrixType=Laplace2D --nx=1000 --ny=1000 --its=20 --scalarCoefficient
    POSTFIX_AND_ARGS_1 "ConstantTensor" --linAlgebra=Tpetra --xml=material.xml --matrixType=Laplace2D --nx=1000 --ny=1000 --its=20 --tensorCoefficient
    POSTFIX_AND_ARGS_2 "TensorFromFile" --linAlgebra=Tpetra --xml=material.xml --matrixType=Laplace2D --nx=100 --ny=100 --its=20 --material=material.mm
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "BlockDiagonalization"
    ARGS "--linAlgebra=Tpetra --xml=block_diagonal_drop.xml --matrixType=Laplace2D --nx=100  --ny=100  --its=20 --blocknumber=blocknumber.mm"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  TRIBITS_ADD_TEST(
    Driver
    NAME "DriverTpetra_Rerun"
    POSTFIX_AND_ARGS_0 "InputDeck" --linAlgebra=Tpetra --xml=scaling-with-rerun.xml
    POSTFIX_AND_ARGS_1 "Repeat" --linAlgebra=Tpetra --reruns=3 --stacked-timer
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

ENDIF()

IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)

  IF(${PACKAGE_NAME}_ENABLE_Zoltan)
    IF(${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra"
        ARGS "--linAlgebra=Tpetra --xml=scaling-complex.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_WithGlobalConstants"
        ARGS "--linAlgebra=Tpetra --xml=scaling-complex-withglobalconstants.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ELSE()
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra"
        ARGS "--linAlgebra=Tpetra --xml=scaling.xml --its=14"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetraSingleReduceCG"
        ARGS "--linAlgebra=Tpetra --xml=scaling.xml --belosType=\"TPETRA CG SINGLE REDUCE\" --its=14"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetraYaml"
        ARGS "--linAlgebra=Tpetra --yaml=scaling.yaml --its=14"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_WithGlobalConstants"
        ARGS "--linAlgebra=Tpetra --xml=scaling-withglobalconstants.xml --its=14"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverEnforceConstraints"
        ARGS "--xml=sa_enforce_constraints.xml --stretchx=3. --matrixType=Laplace2D --nx=43 --ny=40 --its=12"
        NUM_MPI_PROCS 1
        COMM mpi # HAVE_MPI required
        RUN_SERIAL
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverDiagonalModifications"
        ARGS "--matrix=A5submatrix.mm --rhs=rhs5sub.mm --xml=spe10.xml --its=17"
        NUM_MPI_PROCS 1
        COMM mpi # HAVE_MPI required
        RUN_SERIAL
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

      TRIBITS_COPY_FILES_TO_BINARY_DIR(DriverDiagonalModifications_cp SOURCE_FILES A5submatrix.mm rhs5sub.mm spe10.xml)

      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_Distance2Coloring"
        ARGS "--linAlgebra=Tpetra --xml=scaling_distance2_agg.xml --its=15"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )

    ENDIF() # INST_COMPLEX_INT_INT

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetraTripleMatrixProduct"
      ARGS "--linAlgebra=Tpetra --xml=tripleMatrixProduct.xml --its=13"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      PASS_REGULAR_EXPRESSION "Belos converged"
    )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetraILU"
      ARGS "--linAlgebra=Tpetra --xml=sa_with_ilu.xml --its=15"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      PASS_REGULAR_EXPRESSION "Belos converged"
    )

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverTpetraIfpack2LinePartitioner"
      ARGS "--linAlgebra=Tpetra --xml=sa_with_Ifpack2_line_detection.xml --its=48"
      NUM_MPI_PROCS 4
      COMM mpi # HAVE_MPI required
      PASS_REGULAR_EXPRESSION "Belos converged"
    )

  ENDIF() # ENABLE_ZOLTAN

  IF(${PACKAGE_NAME}_ENABLE_Zoltan2Core)
    IF(${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_Milestone"
        ARGS "--linAlgebra=Tpetra --xml=conchas_milestone_zoltan2_complex.xml"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
      )
    ELSE()
      TRIBITS_ADD_TEST(
        Driver
        NAME "DriverTpetra_Milestone"
        ARGS "--linAlgebra=Tpetra --xml=conchas_milestone_zoltan2.xml --its=19"
        NUM_MPI_PROCS 4
        COMM mpi # HAVE_MPI required
        PASS_REGULAR_EXPRESSION "Belos converged"
      )
    ENDIF()

    TRIBITS_ADD_TEST(
      Driver
      NAME "DriverRebalanceInitialSystem"
      ARGS "--linAlgebra=Tpetra --belosType=\"Pseudoblock GMRES\" --tol=1e-8 --xml=scaling.xml --matrix=aux2.mat --coords=coords2.mat --repartXML=repart.xml --its=9"
      NUM_MPI_PROCS 3
      COMM mpi # HAVE_MPI required
      PASS_REGULAR_EXPRESSION "Belos converged"
    )

  ENDIF()

  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "RAPScalingTestTpetra"
    ARGS "--linAlgebra=Tpetra --xml=rap.xml --solver=none --nx=50 --ny=50 --rebuild=1"
    NUM_MPI_PROCS 4
    COMM serial mpi
    #do not check for convergence
  )
  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "CalcRotations"
    ARGS "--linAlgebra=Tpetra --xml=comp_rotations.xml  --matrixType=Elasticity3D --nx=40 --ny=40 --nz=4 --its=26 --belosType=gmres --tol=1.0e-8  --muelu-computed-nullspace"
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Belos converged"
  )
  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "AnisotropicElasticityEmin"
    ARGS "--linAlgebra=Tpetra --xml=emin_with_dl_elasticity3D.xml  --matrixType=Elasticity3D --nx=20 --ny=20 --nz=20 --stretchx=3. --its=30 --belosType=cg --tol=1.0e-8  --muelu-computed-nullspace"
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Belos converged"
  )
  MUELU_ADD_SERIAL_AND_MPI_TEST(
    Driver
    NAME "MaterialAnisotropicElasticityEmin"
    POSTFIX_AND_ARGS_0 "ConstantScalar" --linAlgebra=Tpetra --xml=emin_with_dl_material_elasticity3D.xml --matrixType=Elasticity3D --nx=20 --ny=20 --nz=20 --stretchx=3. --its=30 --belosType=cg
                       --tol=1.0e-8 --muelu-computed-nullspace --scalarCoefficient
    POSTFIX_AND_ARGS_1 "ConstantTensor" --linAlgebra=Tpetra --xml=emin_with_dl_material_elasticity3D.xml --matrixType=Elasticity3D --nx=20 --ny=20 --nz=20 --stretchx=3. --its=30 --belosType=cg
                       --tol=1.0e-8 --muelu-computed-nullspace --tensorCoefficient
    NUM_MPI_PROCS 4
    COMM serial mpi
    PASS_REGULAR_EXPRESSION "Belos converged"
  )

  IF(${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Driver
      NAME "SmootherScalingTestTpetra"
      ARGS "--linAlgebra=Tpetra --xml=smoother_complex.xml --nx=10 --ny=10 --solver=none"
      NUM_MPI_PROCS 4
      COMM serial mpi
      #do not check for convergence
    )
  ELSE()
    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Driver
      NAME "SmootherScalingTestTpetra"
      ARGS "--linAlgebra=Tpetra --xml=smoother.xml --nx=10 --ny=10 --solver=none"
      NUM_MPI_PROCS 4
      COMM serial mpi
      #do not check for convergence
    )
  ENDIF()

  # Reading a matrix saved in real format as complex does not work.
  IF(NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    MUELU_ADD_SERIAL_AND_MPI_TEST(
      ReadMatrix
      NAME "ReadMatrixTpetra"
      ARGS "--linAlgebra=Tpetra --xml=simple.xml --its=10"
      NUM_MPI_PROCS 4
      COMM serial mpi
      PASS_REGULAR_EXPRESSION "Belos converged"
    )
    TRIBITS_ADD_TEST(
      ReadMatrix
      NAME "SmooVecCoalesce"
      ARGS "--linAlgebra=Tpetra --matrixfile=smooVec.mm --xml=smooVecCoalesce.xml --numpdes=3 --its=2"
      NUM_MPI_PROCS 3
      COMM mpi
      PASS_REGULAR_EXPRESSION "Belos converged"
    )
    IF(MUELU_RECURMG)
      MUELU_ADD_SERIAL_AND_MPI_TEST(
        Driver
        NAME "RecursiveMGWithSemiCoarsening"
        ARGS "--linAlgebra=Tpetra --xml=recurMG.xml --matrixType=Elasticity3D --nx=12  --ny=12 --nz=99 --belosType=\"Pseudoblock GMRES\" --noscale --its=50 --tol=1e-5"
        NUM_MPI_PROCS 4
        COMM serial mpi
      )
    ENDIF()

  ENDIF()

ENDIF()

TRIBITS_ADD_EXECUTABLE(
  SpMVPerformance
  SOURCES SpMVPerformance.cpp
  COMM mpi
  CATEGORIES BASIC PERFORMANCE
)

IF(Tpetra_INST_SERIAL)
  TRIBITS_ADD_TEST(
    SpMVPerformance
    NAME_POSTFIX "Tpetra"
    COMM mpi
    NUM_MPI_PROCS 1
    ARGS "--node=serial --linAlgebra=Tpetra --nx=120 --ny=120 --nz=120 --matrixType=Laplace3D --num-runs=1000"
    PASS_REGULAR_EXPRESSION "Complete"
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    SpMVPerformance
    NAME_POSTFIX "Tpetra"
    COMM mpi
    NUM_MPI_PROCS 4
    ARGS "--node=serial --linAlgebra=Tpetra --nx=120 --ny=120 --nz=120 --matrixType=Laplace3D --num-runs=1000"
    PASS_REGULAR_EXPRESSION "Complete"
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )
ENDIF()

IF(${PACKAGE_NAME}_ENABLE_Epetra)
  TRIBITS_ADD_TEST(
    SpMVPerformance
    NAME_POSTFIX "Epetra"
    COMM mpi
    NUM_MPI_PROCS 1
    ARGS "--linAlgebra=Epetra --nx=120 --ny=120 --nz=120 --matrixType=Laplace3D --num-runs=1000"
    PASS_REGULAR_EXPRESSION "Complete"
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )

  TRIBITS_ADD_TEST(
    SpMVPerformance
    NAME_POSTFIX "Epetra"
    COMM mpi
    NUM_MPI_PROCS 4
    ARGS "--linAlgebra=Epetra --nx=120 --ny=120 --nz=120 --matrixType=Laplace3D --num-runs=1000"
    PASS_REGULAR_EXPRESSION "Complete"
    RUN_SERIAL
    CATEGORIES PERFORMANCE
  )
ENDIF()

# Driver to evaluate low level Matrix Matrix multiply kernels potentially using different TPLs
# This differs from the existing Matrix Matrix driver that forms parallel random matrices for testing
# This particular driver does not require MPI and is focused on evaluations using smaller
# MatrixMarket matrices.
IF(${PACKAGE_NAME}_ENABLE_Experimental)

  # Reading a matrix saved in real format as complex does not work.
  IF(NOT ${PACKAGE_NAME}_INST_COMPLEX_INT_INT)
    TRIBITS_ADD_EXECUTABLE(
      MMKernelDriver
      SOURCES MMKernelDriver.cpp
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      MMKernelDriver
      NAME "MatrixMatrixKernelDriver"
      ARGS "--linAlgebra=Tpetra"
      NUM_MPI_PROCS 1
      COMM serial mpi
    )

    TRIBITS_ADD_EXECUTABLE(
      JacobiKernelDriver
      SOURCES JacobiKernelDriver.cpp
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      JacobiKernelDriver
      NAME "JacobiKernelDriver"
      ARGS "--linAlgebra=Tpetra"
      NUM_MPI_PROCS 1
      COMM serial mpi
    )

    TRIBITS_ADD_EXECUTABLE(
      TwoMatrixMMKernelDriver
      SOURCES TwoMatrixMMKernelDriver.cpp
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      TwoMatrixMMKernelDriver
      NAME "TwoMatrixMMKernelDriver"
      ARGS "--linAlgebra=Tpetra"
      NUM_MPI_PROCS 1
      COMM serial mpi
    )

  ENDIF()
ENDIF()

ASSERT_DEFINED(
  ${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT
  Tpetra_INST_SERIAL
  Tpetra_INST_CUDA
  Tpetra_INST_HIP
)

IF(${PACKAGE_NAME}_INST_DOUBLE_INT_LONGLONGINT)

  TRIBITS_ADD_EXECUTABLE(
    MatvecKernelDriver
    SOURCES MatvecKernelDriver.cpp
    COMM serial mpi
    CATEGORIES BASIC PERFORMANCE
  )

  IF(Tpetra_INST_SERIAL)

    SET(COUNTER 1)
    SET(MAX_COUNT 11)

    WHILE(COUNTER LESS MAX_COUNT)
      MATH(EXPR n "150 * ${COUNTER}")
      MATH(EXPR NP "${COUNTER} * ${COUNTER}")

      TRIBITS_ADD_TEST(
        MatvecKernelDriver
        NAME Matvec_Serial
        COMM mpi
        ARGS "--stackedtimer --matrixType=Laplace2D --nx=${n} --ny=${n} --node=serial"
        NUM_MPI_PROCS ${NP}
        RUN_SERIAL
        CATEGORIES PERFORMANCE
      )

      MATH(EXPR COUNTER "${COUNTER} + 1")
    ENDWHILE()

  ENDIF()

  IF(Tpetra_INST_CUDA)

    SET(COUNTER 1)
    SET(MAX_COUNT 11)

    WHILE(COUNTER LESS MAX_COUNT)
      MATH(EXPR n "600 * ${COUNTER}")
      MATH(EXPR NP "${COUNTER} * ${COUNTER}")

      TRIBITS_ADD_TEST(
        MatvecKernelDriver
        NAME Matvec_CUDA
        COMM mpi
        ARGS "--stackedtimer --matrixType=Laplace2D --nx=${n} --ny=${n} --node=cuda"
        NUM_MPI_PROCS ${NP}
        RUN_SERIAL
        CATEGORIES PERFORMANCE
      )

      MATH(EXPR COUNTER "${COUNTER} + 1")
    ENDWHILE()

  ENDIF()

  IF(Tpetra_INST_HIP)

    SET(COUNTER 1)
    SET(MAX_COUNT 11)

    WHILE(COUNTER LESS MAX_COUNT)
      MATH(EXPR n "600 * ${COUNTER}")
      MATH(EXPR NP "${COUNTER} * ${COUNTER}")

      TRIBITS_ADD_TEST(
        MatvecKernelDriver
        NAME Matvec_HIP
        COMM mpi
        ARGS "--stackedtimer --matrixType=Laplace2D --nx=${n} --ny=${n} --node=hip"
        NUM_MPI_PROCS ${NP}
        RUN_SERIAL
        CATEGORIES PERFORMANCE
      )

      MATH(EXPR COUNTER "${COUNTER} + 1")
    ENDWHILE()

  ENDIF()

  IF(Tpetra_INST_SYCL)

    SET(COUNTER 1)
    SET(MAX_COUNT 11)

    WHILE(COUNTER LESS MAX_COUNT)
      MATH(EXPR n "600 * ${COUNTER}")
      MATH(EXPR NP "${COUNTER} * ${COUNTER}")

      TRIBITS_ADD_TEST(
        MatvecKernelDriver
        NAME Matvec_SYCL
        COMM mpi
        ARGS "--stackedtimer --matrixType=Laplace2D --nx=${n} --ny=${n} --node=sycl"
        NUM_MPI_PROCS ${NP}
        RUN_SERIAL
        CATEGORIES PERFORMANCE
      )

      MATH(EXPR COUNTER "${COUNTER} + 1")
    ENDWHILE()

  ENDIF()
ENDIF()
