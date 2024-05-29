# Ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Belos
  ${PACKAGE_NAME}_ENABLE_ML
  ${PACKAGE_NAME}_ENABLE_Zoltan2Core
  )

TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Transfers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../gallery)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../adapters/belos)
TRIBITS_INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)


IF (${PACKAGE_NAME}_ENABLE_Amesos2)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ReitzingerPFactory
    SOURCES ReitzingerPFactory.cpp
    ARGS ""
    COMM serial mpi
    )
ENDIF()


IF (${PACKAGE_NAME}_ENABLE_Belos)

  TRIBITS_ADD_EXECUTABLE(
    Maxwell3D
    SOURCES Maxwell3D.cpp
    COMM serial mpi
    )

  IF (${PACKAGE_NAME}_ENABLE_Amesos2)
    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-Tpetra"
      ARGS "--linAlgebra=Tpetra --reuse --xml=Maxwell.xml"
      ARGS "--linAlgebra=Tpetra --precType=MueLu-Maxwell1 --xml=Maxwell_Reitzinger.xml"
      ARGS "--linAlgebra=Tpetra --precType=MueLu-Maxwell1 --xml=Maxwell_Reitzinger_repart.xml"
      ARGS "--linAlgebra=Tpetra --precType=MueLu-Maxwell1 --xml=Maxwell_Reitzinger_repart.xml --ensure-kn"
      COMM serial mpi
      NUM_MPI_PROCS 4
      )

    IF (${PACKAGE_NAME}_ENABLE_Zoltan2Core)
      TRIBITS_ADD_TEST(
        Maxwell3D
        NAME "Maxwell3D-Tpetra2"
        ARGS "--linAlgebra=Tpetra --reuse --xml=Maxwell2.xml"
        COMM serial mpi
        NUM_MPI_PROCS 4
      )
    ENDIF()

    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-GMHD"
      ARGS "--linAlgebra=Tpetra --xml=gmhd.xml --SM=GMHD_curlcurlBlk.dat --D0=GMHD_D0.dat --precType=MueLu-Maxwell1  --belos_xml=Belos.xml --tol=1e9 --coords=GMHD_nodalCoords.dat --GmhdA=GMHD.dat --belosSolverType=\"Block GMRES\" "
      COMM serial mpi
      NUM_MPI_PROCS 1
      )

    TRIBITS_ADD_TEST(
      Maxwell3D
      NAME "Maxwell3D-GMHD"
      ARGS "--linAlgebra=Tpetra --xml=gmhd.xml --SM=GMHD_curlcurlBlk.dat --D0=GMHD_D0.dat --precType=MueLu-Maxwell1  --belos_xml=Belos.xml --tol=1e9 --coords=GMHD_nodalCoords.dat --GmhdA=GMHD.dat --belosSolverType=\"Block GMRES\" "
      COMM serial mpi
      NUM_MPI_PROCS 2
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
    SOURCE_FILES M0.mat M1.mat S.mat D0.mat coords.mat Belos.xml Stratimikos.xml Stratimikos_Maxwell1.xml Maxwell.xml Maxwell2.xml Maxwell3.xml Maxwell_ML.xml Maxwell_ML_MueLu.xml Maxwell_ML1.xml Maxwell_ML1_MueLu.xml Maxwell_Reitzinger.xml Maxwell_Reitzinger_repart.xml GMHD.dat GMHD_nodalCoords.dat  GMHD_D0.dat GMHD_curlcurlBlk.dat gmhd.xml
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(Maxwell_cp_complex
    SOURCE_FILES M0_complex.mat M1_complex.mat S_complex.mat D0_complex.mat Maxwell_complex.xml
    )


ENDIF()
