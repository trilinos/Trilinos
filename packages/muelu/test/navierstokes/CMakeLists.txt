INCLUDE(TribitsAddExecutableAndTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Graph) # we need this line because this file is not on the header file List of MueLu
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../../src/Smoothers)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_EpetraExt
  ${PACKAGE_NAME}_ENABLE_Amesos
  )


IF (${PACKAGE_NAME}_ENABLE_Epetra AND
    ${PACKAGE_NAME}_ENABLE_EpetraExt AND
    ${PACKAGE_NAME}_ENABLE_Ifpack AND
    ${PACKAGE_NAME}_ENABLE_Amesos)

  IF(Trilinos_ENABLE_AztecOO)
    # build only if AztecOO available, too
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2D_Epetra
      SOURCES Navier2D_Epetra.cpp
      COMM mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_Epetra
      SOURCES Navier2DBlocked_Epetra.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_xml
      SOURCES Navier2DBlocked_xml.cpp
      ARGS "--xml=myXML.xml"
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_xml2
      SOURCES Navier2DBlocked_xml.cpp
      ARGS "--xml=myXMLg.xml"
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_BraessSarazinSmoother
      SOURCES Navier2DBlocked_BraessSarazinSmoother.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    #TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #  Navier2DBlocked_SimpleSmoother
    #  SOURCES Navier2DBlocked_SimpleSmoother.cpp
    #  COMM serial mpi
    #  DEPLIBS muelu xpetra
    #  )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_BraessSarazin
      SOURCES Navier2DBlocked_BraessSarazin.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    # tests with rebalancing
    TRIBITS_ADD_EXECUTABLE(
      Navier2DBlocked_rebalancing
      SOURCES Navier2DBlocked_test.cpp
      COMM mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_TEST(
      Navier2DBlocked_rebalancing
      NAME "Navier2DBlocked_reb_1_0"
      ARGS "--rebalanceBlock0=1 --rebalanceBlock1=0"
      NUM_MPI_PROCS 4
      COMM mpi
      )
    TRIBITS_ADD_TEST(
      Navier2DBlocked_rebalancing
      NAME "Navier2DBlocked_reb_0_1"
      ARGS "--rebalanceBlock0=0 --rebalanceBlock1=1"
      NUM_MPI_PROCS 4
      COMM mpi
      )
    TRIBITS_ADD_TEST(
      Navier2DBlocked_rebalancing
      NAME "Navier2DBlocked_reb_1_1"
      ARGS "--rebalanceBlock0=1 --rebalanceBlock1=1"
      NUM_MPI_PROCS 4
      COMM mpi
      )

    TRIBITS_ADD_EXECUTABLE(
      Navier2DBlocked_rebalancing2
      SOURCES Navier2DBlocked_test2.cpp
      COMM mpi
      DEPLIBS muelu xpetra
      )
      
    TRIBITS_ADD_TEST(
      Navier2DBlocked_rebalancing2
      NAME "Navier2DBlockedReuseAggs"
      ARGS "--rebalanceBlocks=0"
      NUM_MPI_PROCS 4
      COMM mpi
      )
    TRIBITS_ADD_TEST(
      Navier2DBlocked_rebalancing2
      NAME "Navier2DBlockedReuseAggsReb"
      ARGS "--rebalanceBlocks=1"
      NUM_MPI_PROCS 4
      COMM mpi
      )
            
    #TRIBITS_ADD_EXECUTABLE_AND_TEST(
    #  Navier2DBlocked_Simple
    #  SOURCES Navier2DBlocked_Simple.cpp
    #  COMM serial mpi
    #  DEPLIBS muelu xpetra
    #  )

    TRIBITS_COPY_FILES_TO_BINARY_DIR(Navier2D_cp
      SOURCE_FILES A5932_re1000.txt b5932_re1000.txt myXML.xml myXMLg.xml
      )
  ENDIF()

ENDIF()
