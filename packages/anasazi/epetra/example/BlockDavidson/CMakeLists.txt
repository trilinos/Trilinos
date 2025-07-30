

IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    Epetra_BD_Galeri_example
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




