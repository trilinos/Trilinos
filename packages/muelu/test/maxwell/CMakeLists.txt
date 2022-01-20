# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_ML
  )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)


 TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ReitzingerPFactory
    SOURCES ReitzingerPFactory.cpp
    ARGS ""
    COMM serial mpi
    )



IF (${PACKAGE_NAME}_ENABLE_Belos AND ${PACKAGE_NAME}_ENABLE_Kokkos_Refactor)

  TRIBITS_ADD_EXECUTABLE(
    Maxwell3D
    SOURCES Maxwell3D.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Tpetra)

    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-Tpetra"
      ARGS "--linAlgebra=Tpetra --reuse --xml=Maxwell.xml"
      ARGS "--linAlgebra=Tpetra --reuse --xml=Maxwell2.xml"
      ARGS "--linAlgebra=Tpetra --precType=MueLu-Maxwell1 --xml=Maxwell_Reitzinger.xml"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    IF (${PACKAGE_NAME}_ENABLE_ML AND ${PACKAGE_NAME}_ENABLE_Epetra)
      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Tpetra-ML-RefMaxwell-list"
        ARGS "--linAlgebra=Tpetra --xml=Maxwell_ML_MueLu.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Tpetra-ML-Maxwell1-list"
        ARGS "--linAlgebra=Tpetra --precType=MueLu-Maxwell1 --xml=Maxwell_ML1_MueLu.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

    ENDIF()

    IF (${PACKAGE_NAME}_INST_DOUBLE_INT_INT AND ${PACKAGE_NAME}_ENABLE_Stratimikos)
      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Tpetra-Stratimikos"
        ARGS "--linAlgebra=Tpetra --solverName=Stratimikos --instantiation=DOUBLE_INT_INT --reuse --xml=Maxwell.xml"
        ARGS "--linAlgebra=Tpetra --solverName=Stratimikos --instantiation=DOUBLE_INT_INT --reuse --xml=Maxwell2.xml"
        ARGS "--linAlgebra=Tpetra --solverName=Stratimikos --instantiation=DOUBLE_INT_INT --reuse --xml=Maxwell3.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )
    ENDIF()
  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Epetra)

    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-Epetra"
      ARGS "--linAlgebra=Epetra --reuse --xml=Maxwell.xml"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    IF (${PACKAGE_NAME}_ENABLE_Stratimikos)
      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Epetra-Stratimikos"
        ARGS "--linAlgebra=Epetra --solverName=Stratimikos --reuse --xml=Maxwell.xml"
        ARGS "--linAlgebra=Epetra --solverName=Stratimikos --reuse --xml=Maxwell3.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )
    ENDIF()

    IF (${PACKAGE_NAME}_ENABLE_ML)

      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-ML-RefMaxwell"
        ARGS "--linAlgebra=Epetra --precType=ML-RefMaxwell  --xml=Maxwell_ML.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-ML-Maxwell1"
        ARGS "--linAlgebra=Epetra --precType=ML-Maxwell  --xml=Maxwell_ML1.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
        )

      # Hiptmair not implemented for MueLu+Epetra
      # TRIBITS_ADD_TEST(
      #   Maxwell3D
      #   NAME "Maxwell3D-ML-list"
      #   ARGS "--linAlgebra=Epetra --precType=MueLu-Maxwell1  --xml=Maxwell_ML1_muelu.xml"
      #   COMM serial mpi
      #   NUM_MPI_PROCS 4
      #   )
    ENDIF()
  ENDIF()


  TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp
    SOURCE_FILES M0.mat M1.mat S.mat D0.mat coords.mat Belos.xml Stratimikos.xml Maxwell.xml Maxwell2.xml Maxwell3.xml Maxwell_ML.xml Maxwell_ML_MueLu.xml Maxwell_ML1.xml Maxwell_ML1_MueLu.xml Maxwell_Reitzinger.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp_complex
    SOURCE_FILES M0_complex.mat M1_complex.mat S_complex.mat D0_complex.mat Maxwell_complex.xml
    )


ENDIF()
