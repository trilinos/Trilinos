INCLUDE(XpetraAddTest)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

IF (${PACKAGE_NAME}_ENABLE_Galeri)

#
# Stratimikos example
#

IF((${PACKAGE_NAME}_ENABLE_Stratimikos) AND
   (${PACKAGE_NAME}_ENABLE_Epetra AND ${PACKAGE_NAME}_ENABLE_EpetraExt AND
    ${PACKAGE_NAME}_ENABLE_Ifpack AND ${PACKAGE_NAME}_ENABLE_Amesos))

  TRIBITS_ADD_EXECUTABLE(
    Stratimikos
    SOURCES Stratimikos.cpp
    DEPLIBS muelu xpetra muelu-adapters
    COMM serial mpi
    )

  TRIBITS_COPY_FILES_TO_BINARY_DIR(advanced_stratimikos_cp
    SOURCE_FILES stratimikos_ParameterList.xml stratimikos_ParameterList2.xml
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Stratimikos
    NAME Stratimikos2
    ARGS "--xml=stratimikos_ParameterList2.xml"
    NUM_MPI_PROCS 4
    COMM serial mpi
    )

ENDIF()

ENDIF()
