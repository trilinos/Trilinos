INCLUDE(PackageAddExecutableAndTest)

PACKAGE_ADD_EXECUTABLE(
  Epetra_BD_example
  SOURCES BlockDavidsonEpetraEx.cpp
  COMM serial mpi
)

PACKAGE_ADD_EXECUTABLE(
  Epetra_BD_Gen_example
  SOURCES BlockDavidsonEpetraExGen.cpp
  COMM serial mpi
)

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BD_File_example
    SOURCES BlockDavidsonEpetraExFile.cpp
    COMM serial mpi
  )

ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)


IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  PACKAGE_ADD_EXECUTABLE(
    Epetra_BD_Ifpack_example
    SOURCES BlockDavidsonEpetraExGenPrecIfpack.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
