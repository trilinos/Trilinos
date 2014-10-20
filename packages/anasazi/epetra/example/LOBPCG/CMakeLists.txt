
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

ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_EpetraExt)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Triutils)
ASSERT_DEFINED(${PACKAGE_NAME}_ENABLE_Ifpack)

IF (${PACKAGE_NAME}_ENABLE_EpetraExt)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_LOBPCG_shifted_example
    SOURCES LOBPCGEpetraExGenShifted.cpp
    COMM serial mpi
  )

IF (${PACKAGE_NAME}_ENABLE_Triutils)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_LOBPCG_File_example
    SOURCES LOBPCGEpetraExFile.cpp
    COMM serial mpi
  )

ENDIF(${PACKAGE_NAME}_ENABLE_Triutils)
ENDIF(${PACKAGE_NAME}_ENABLE_EpetraExt)


IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  TRIBITS_ADD_EXECUTABLE(
    Epetra_LOBPCG_Ifpack_example
    SOURCES LOBPCGEpetraExGenPrecIfpack.cpp
    COMM serial mpi
    )

ENDIF(${PACKAGE_NAME}_ENABLE_Ifpack)
