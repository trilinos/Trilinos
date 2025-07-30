
IF (${PACKAGE_NAME}_ENABLE_Galeri)
  TRIBITS_ADD_EXECUTABLE(
    Epetra_LOBPCG_Galeri_example
    SOURCES LOBPCGEpetra.cpp
    COMM serial mpi
  )
ENDIF ()

TRIBITS_ADD_EXECUTABLE(
  Epetra_LOBPCG_example
  SOURCES LOBPCGEpetraEx.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Epetra_LOBPCG_simple_example
  SOURCES LOBPCGEpetraExSimple.cpp
  COMM serial mpi
)

TRIBITS_ADD_EXECUTABLE(
  Epetra_LOBPCG_Gen_example
  SOURCES LOBPCGEpetraExGen.cpp
  COMM serial mpi
)




