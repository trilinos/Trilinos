

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    BlockDavidsonEpetra
    SOURCES BlockDavidsonEpetra.cpp
    COMM serial mpi
  )
ENDIF ()


TRIBITS_ADD_EXECUTABLE(
  Epetra_BD_example
  SOURCES BlockDavidsonEpetraEx.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Epetra_BD_Gen_example
  SOURCES BlockDavidsonEpetraExGen.cpp
  COMM serial mpi
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_BD_Shift_example
    SOURCES BlockDavidsonEpetraExGenShifted.cpp
    COMM serial mpi
  )

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_BD_File_example
    SOURCES BlockDavidsonEpetraExFile.cpp
    COMM serial mpi
  )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)


IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_BD_Ifpack_example
    SOURCES BlockDavidsonEpetraExGenPrecIfpack.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
