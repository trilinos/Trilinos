
INCLUDE(Trilinos_Add_Executable_And_Test)
INCLUDE(Copy_Files_To_Binary_Dir)


#
# simple_stratimikos_example
#


TRILINOS_ADD_EXECUTABLE(
  simple_stratimikos_example
  SOURCES
    simple_stratimikos_example.cpp
  COMM serial mpi
  )


COPY_FILES_TO_BINARY_DIR(CopyExampleFourByFourMtx
  DEST_FILES FourByFour.mtx
  SOURCE_DIR ${PACKAGE_SOURCE_DIR}/../amesos/test/Test_Basic
  SOURCE_PREFIX "In_"
  EXEDEPS simple_stratimikos_example
  )


COPY_FILES_TO_BINARY_DIR(CopyExampleXMLFiles
  DEST_FILES
    amesos.klu.xml
    aztecoo.ifpack.xml
    aztecoo.ml.xml
    belos.ifpack.xml
    belos.ml.xml
  SOURCE_PREFIX "_"
  EXEDEPS simple_stratimikos_example
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Amesos)
IF (${PACKAGE_NAME}_ENABLE_Amesos)
  TRILINOS_ADD_TEST(
    simple_stratimikos_example
    NAME simple_stratimikos_example_amesos
    ARGS 
      "--linear-solver-params-file=amesos.klu.xml --matrix-file=FourByFour.mtx --tol=1e-13"
      "--extra-linear-solver-params=\"<ParameterList><Parameter name='Linear Solver Type' type='string' value='Amesos'/></ParameterList>\" --matrix-file=FourByFour.mtx --tol=1e-10 --linear-solver-params-used-file=lowsf.amesos.used.xml"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )
ENDIF()


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
IF (${PACKAGE_NAME}_ENABLE_AztecOO)


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
  IF (${PACKAGE_NAME}_ENABLE_Ifpack)
    TRILINOS_ADD_TEST(
      simple_stratimikos_example
      NAME simple_stratimikos_example_ifpack
      ARGS 
        "--linear-solver-params-file=aztecoo.ifpack.xml --matrix-file=FourByFour.mtx --tol=1e-13"
        "--extra-linear-solver-params=\"<ParameterList><Parameter name='Linear Solver Type' type='string' value='AztecOO'/></ParameterList>\" --matrix-file=FourByFour.mtx --tol=1e-10 --linear-solver-params-used-file=lowsf.aztecoo.used.xml"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  ENDIF()


  ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_ML)
  IF (${PACKAGE_NAME}_ENABLE_ML)
    TRILINOS_ADD_TEST(
      simple_stratimikos_example
      NAME simple_stratimikos_example_ml
      ARGS 
        "--linear-solver-params-file=aztecoo.ml.xml --matrix-file=FourByFour.mtx --tol=1e-13"
      COMM serial mpi
      NUM_MPI_PROCS 1
      )
  ENDIF()


ENDIF()

# ToDo: Add Belos tests!


#
# Generate sample output
#


# ToDo: Implement!


#
# MixedOrderPhysicsBasedPreconditioner
#


TRILINOS_ADD_EXECUTABLE(
  MixedOrderPhysicsBasedPreconditioner
  SOURCES
    MixedOrderPhysicsBasedPreconditioner.cpp
  COMM serial mpi
  )


COPY_FILES_TO_BINARY_DIR(CopyExampleMOPBPFiles
  DEST_FILES
    MixedOrderPhysicsBasedPreconditioner.Base.xml
    MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml
    MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml
    MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml
    MixedOrderPhysicsBasedPreconditioner.Belos.xml
    P2.mtx
    P1.mtx
    M11.mtx
    M12.mtx
    M21.mtx
    M22.mtx
  SOURCE_PREFIX "_"
  EXEDEPS MixedOrderPhysicsBasedPreconditioner
  )


ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_AztecOO)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)
IF (${PACKAGE_NAME}_ENABLE_AztecOO AND ${PACKAGE_NAME}_ENABLE_Ifpack)

  TRILINOS_ADD_TEST(
    MixedOrderPhysicsBasedPreconditioner
    NAME MixedOrderPhysicsBasedPreconditioner_aztecoo_ifpack
    ARGS 
      "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --solve-tol=1.0"
      "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml --solve-tol=1.0"
      "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.AllSolveOutput.xml --solve-tol=1.0"
      "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutputMSolveOutputToFile.xml --solve-tol=1.0"
      "--params-file=MixedOrderPhysicsBasedPreconditioner.Base.xml --extra-params-file=MixedOrderPhysicsBasedPreconditioner.P2SolveOutput.xml --use-algebraic-prec --solve-tol=1.0"
    COMM serial mpi
    NUM_MPI_PROCS 1
    )

ENDIF()


# ToDo: Add Belos tests!
