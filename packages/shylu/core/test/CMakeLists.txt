
#TRIBITS_ADD_EXECUTABLE_AND_TEST(
#    aztec_driver
#    SOURCES shylu_driver.cpp
#    COMM mpi
#    )

TRIBITS_ADD_EXECUTABLE_AND_TEST(
    belos_driver
    SOURCES shylu_belos_driver.cpp
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

TRIBITS_ADD_EXECUTABLE(
    iqr_driver
    SOURCES shylu_iqr_driver.cpp
    COMM serial mpi
    )

TRIBITS_ADD_EXECUTABLE(
	epetra_interface
	SOURCES epetra_interface_test.cpp 
	COMM serial mpi
	)

IF(HAVE_SHYLUCORE_TPETRA)

  TRIBITS_ADD_EXECUTABLE(
	tpetra_interface
	SOURCES tpetra_interface_test.cpp 
	COMM serial mpi
	)

ENDIF()

IF(HAVE_SHYLUCORE_TPETRA AND HAVE_SHYLUCORE_AMESOS2 AND TPL_ENABLE_SuperLU AND ${PACKAGE_NAME}_ENABLE_ParMETIS)
  TRIBITS_ADD_TEST(
	tpetra_interface
        NAME tpetra_interface_with_superlu
	COMM serial mpi
	NUM_MPI_PROCS 4
	)
ENDIF()

IF (${PACKAGE_NAME}_ENABLE_ParMETIS) 
  TRIBITS_ADD_TEST(
    epetra_interface
    NAME epetra_interface_with_parmetis
    COMM serial mpi
    NUM_MPI_PROCS 4
    )
ENDIF ()


IF (${PACKAGE_NAME}_ENABLE_ParMETIS)
   TRIBITS_ADD_TEST(
   iqr_driver
   NAME iqr_driver_with_parmetis
   COMM serial mpi
   NUM_MPI_PROCS 4
    )
ENDIF()

TRIBITS_COPY_FILES_TO_BINARY_DIR(ShyLUDriverFiles
    SOURCE_FILES ShyLU.xml shylu_iqr_parameters.xml
    shylu_no_iqr_parameters.xml wathen96.mtx wathenSmall.mtx wathenTiny.mtx
    shylu_iqr_parameters.seq.xml
    ShyLU_epetra_interface.xml
    EXEDEPS belos_driver iqr_driver epetra_interface
    )

