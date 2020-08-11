INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/unit_tests)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../test/scaling)

#
# Executables
#

IF (${PACKAGE_NAME}_ENABLE_Tpetra AND ${PACKAGE_NAME}_ENABLE_Ifpack2 AND ${PACKAGE_NAME}_ENABLE_Amesos2)
  SET(${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK YES)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos)
  SET(${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK YES)
ENDIF()



IF (${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)
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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Simple_cp
    SOURCE_FILES simple.xml simple.yaml
    )

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

  TRIBITS_COPY_FILES_TO_BINARY_DIR(advanced_stratimikos_cp
    SOURCE_FILES stratimikos_ParameterList.xml stratimikos_ParameterList.yaml stratimikos_ParameterList2.xml
    set1-gmres.xml set1-cg.xml
    set2-sgs1.xml set2-sgs3.xml
    set3-mg-jacobi.xml set3-mg-sgs.xml set3-mg-gs.xml set3-mg-chebyshev.xml
    mg-gpu.xml
    )

  IF (${PACKAGE_NAME}_HAVE_EPETRA_SOLVER_STACK)

    TRIBITS_ADD_TEST(
      Stratimikos
      NAME Stratimikos
      ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Epetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

  ENDIF()

  IF (${PACKAGE_NAME}_HAVE_TPETRA_SOLVER_STACK)

    TRIBITS_ADD_TEST(
      Stratimikos
      NAME Stratimikos
      ARGS "--xml=stratimikos_ParameterList2.xml --linAlgebra=Tpetra"
      NUM_MPI_PROCS 4
      COMM serial mpi
      )

  ENDIF()

ENDIF()
