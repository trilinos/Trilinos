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
    ${PACKAGE_NAME}_ENABLE_Amesos) 

  IF(Trilinos_ENABLE_AztecOO) 
    # build only if AztecOO available, too
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2D_Test
      SOURCES Navier2D_epetra.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )
    
    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_Test
      SOURCES Navier2D_blocked_epetra.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_BraessSarazinSmoother
      SOURCES Navier2D_blocked_bssmoother.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )

    TRIBITS_ADD_EXECUTABLE_AND_TEST(
      Navier2DBlocked_BraessSarazin
      SOURCES Navier2D_blocked_bs.cpp
      COMM serial mpi
      DEPLIBS muelu xpetra
      )
    TRIBITS_COPY_FILES_TO_BINARY_DIR(Navier2D_test_copyfiles
      SOURCE_FILES A5932_re1000.txt b5932_re1000.txt
      )
  ENDIF()
  
ENDIF()
