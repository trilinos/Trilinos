INCLUDE(XpetraAddTest)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/scaling)

#
# Executables
#

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

IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)
  TRIBITS_ADD_EXECUTABLE(
    Simple
    SOURCES Simple.cpp
    COMM serial mpi
  )

  TRIBITS_ADD_TEST(
    Simple
    NAME "SimpleTpetra"
    ARGS "--linAlgebra=Tpetra --xml=simple.xml"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

  TRIBITS_ADD_TEST(
    Simple
    NAME "SimpleTpetraYaml"
    ARGS "--linAlgebra=Tpetra --yaml=simple.yaml"
    NUM_MPI_PROCS 4
    COMM mpi # HAVE_MPI required
  )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple_cp SOURCE_FILES simple.xml simple.yaml)

ENDIF()

#
# Stratimikos example
#
IF(${PACKAGE_NAME}_ENABLE_Stratimikos AND ${PACKAGE_NAME}_ENABLE_Thyra)

  # build executable
  TRIBITS_ADD_EXECUTABLE(
    Stratimikos
    SOURCES Stratimikos.cpp
    COMM serial mpi
  )

  # Teko build
  IF(${PACKAGE_NAME}_ENABLE_Teko AND ${PACKAGE_NAME}_ENABLE_Ifpack2)
    TRIBITS_ADD_EXECUTABLE(
      Teko
      SOURCES Teko.cpp
      COMM serial mpi
    )
  ENDIF()

  TRIBITS_COPY_FILES_TO_BINARY_DIR(
    advanced_stratimikos_cp
    SOURCE_FILES
      mg-gpu.xml
      set1-cg.xml
      set1-gmres.xml
      set2-sgs1.xml
      set2-sgs3.xml
      set3-mg-chebyshev.xml
      set3-mg-gs.xml
      set3-mg-jacobi.xml
      set3-mg-sgs-reuse.xml
      set3-mg-sgs.xml
      stratimikos_ParameterList.xml
      stratimikos_ParameterList.yaml
      stratimikos_ParameterList2.xml
  )

  IF(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Stratimikos
      NAME Stratimikos
      ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Stratimikos
      NAME HandsOn
      POSTFIX_AND_ARGS_0 "GMRES" "--xml=set1-gmres.xml"
      POSTFIX_AND_ARGS_1 "CG" "--xml=set1-cg.xml"
      WILL_FAIL
      NUM_MPI_PROCS 1
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Stratimikos
      NAME HandsOn
      POSTFIX_AND_ARGS_0 "SGS1" "--xml=set2-sgs1.xml"
      POSTFIX_AND_ARGS_1 "SGS3" "--xml=set2-sgs3.xml"
      POSTFIX_AND_ARGS_2 "MG-Jacobi" "--xml=set3-mg-jacobi.xml"
      POSTFIX_AND_ARGS_3 "MG-SGS" "--xml=set3-mg-sgs.xml"
      POSTFIX_AND_ARGS_4 "MG-GS" "--xml=set3-mg-gs.xml"
      POSTFIX_AND_ARGS_5 "MG-Chebyshev" "--xml=set3-mg-chebyshev.xml"
      POSTFIX_AND_ARGS_6 "MG-SGS-reuse" "--xml=set3-mg-sgs-reuse.xml"
      NUM_MPI_PROCS 1
      COMM serial mpi
    )

    MUELU_ADD_SERIAL_AND_MPI_TEST(
      Stratimikos
      NAME Stratimikos_direct
      ARGS "--xml=set3-mg-jacobi.xml --linAlgebra=Tpetra --useDirectMatrixConstruction --matrixType=Laplace1D --nx=3000"
      NUM_MPI_PROCS 4
      COMM serial mpi
    )

  ENDIF()

ENDIF()
